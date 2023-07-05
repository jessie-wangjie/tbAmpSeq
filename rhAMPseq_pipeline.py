#!/usr/bin/env python
"""
Tome rhAMP-seq analysis
"""

import argparse 
import glob
import datetime
from utils.base import *
from utils.common_functions import *
import subprocess
import sys
import pandas as pd
import numpy as np
import io

def window_quantification(cs2_folder, quantification_windows):
    # Amplicon:Window_name:Window_region:flanking_bp. 
    # Bp positions in the amplicon sequence specifying the quantification window, 1-index
    try:
        cs2_info = CRISPRessoShared.load_crispresso_info(cs2_folder)
    except Exception:
        return {}

    if not cs2_info["running_info"]["args"].write_detailed_allele_table:
        raise Exception('CRISPResso run must be run with the parameter --write_detailed_allele_table')

    z = zipfile.ZipFile(os.path.join(cs2_folder, cs2_info["running_info"]["allele_frequency_table_zip_filename"]))
    zf = z.open(cs2_info["running_info"]["allele_frequency_table_filename"])
    df_alleles = pd.read_csv(zf, sep="\t")
    df_alleles["ref_positions"] = df_alleles["ref_positions"].apply(arrstr_to_arr)
    
    # generate the stats JSON for the result schema
    b_json = {"total_read_num": CRISPRessoCORE.get_n_reads_fastq(cs2_info["running_info"]["args"].fastq_r1),
              "merged_r1r2_read_num": int(cs2_info["running_info"]["alignment_stats"]["N_TOT_READS"]),
              "wt_aligned_read_num": int(cs2_info["results"]["alignment_stats"]["counts_total"]['Reference'])}
    if "Beacon" in cs2_info["results"]["alignment_stats"]["counts_total"]:
        b_json["beacon_aligned_read_num"] = int(cs2_info["results"]["alignment_stats"]["counts_total"]["Beacon"])
    elif "PE" in cs2_info["results"]["alignment_stats"]["counts_total"]:
        b_json["beacon_aligned_read_num"] = int(cs2_info["results"]["alignment_stats"]["counts_total"]["PE"])
    else:
        b_json["beacon_aligned_read_num"] = 0

    b_json["aligned_percentage"] = format(100 * (b_json["wt_aligned_read_num"] + b_json["beacon_aligned_read_num"]) / b_json["merged_r1r2_read_num"],
                                          ".1f")
    b_json["wt_aligned_percentage"] = format(
        100 * b_json["wt_aligned_read_num"] / (b_json["wt_aligned_read_num"] + b_json["beacon_aligned_read_num"]), ".1f")
    b_json["beacon_placement_percentage"] = format(
        100 * b_json["beacon_aligned_read_num"] / (b_json["wt_aligned_read_num"] + b_json["beacon_aligned_read_num"]), ".1f")

    qw_stats = []
    for window in quantification_windows:  
        ref_name, qw_name, qw, flank_bp = window.split(":")
        df_alleles["Reference_Name"] = ref_name
        start, end = qw.split("-")

        stats = {"amplicon": ref_name, "window_name": qw_name, "window_region": qw + ":" + flank_bp}
        df_ref = df_alleles[df_alleles["Reference_Name"] == ref_name]
        if df_ref.empty:
            continue
        
        df = df_ref.apply(lambda row: get_modified_in_quantification_window(row, set(range(int(start) - 1, int(end)))), axis=1, result_type='expand')
        g = df.groupby("classification").sum()
        for i in g.index:
            stats[i] = g.loc[i]["#Reads"]
            if i == "modified":
                stats.update(g.loc[i][["indels", "insertion", "deletion", "substitution", "whole_window_deletion"]])

        if int(flank_bp):
            include_idx = []
            include_idx.extend(range(int(start) - int(flank_bp) - 1, int(start) - 1))
            include_idx.extend(range(int(end), int(end) + int(flank_bp)))
            df_flank = pd.concat([df, df_ref.apply(lambda row: get_modified_in_quantification_window(row, sorted(include_idx)), axis=1,
                                                   result_type='expand').add_suffix("_flank").drop("#Reads_flank", axis=1)], axis=1)
            g = df_flank.groupby(["classification", "classification_flank"]).sum()
            for i, j in g.index:
                stats[i + "_" + j + "_flank"] = g.loc[i, j]["#Reads"]
                if j == "modified":
                    stats.update(g.loc[i, j][g.columns.str.endswith("_flank")].add_prefix(i + "_"))

        qw_stats.append(stats)

        if ref_name == "Beacon" and qw_name == "beacon_whole":
            if "indels" in stats:
                b_json["beacon_indel_read_num"] = int(stats["indels"])
                b_json["beacon_sub_read_num"] = int(stats["substitution"])
            else:
                b_json["beacon_indel_read_num"] = 0
                b_json["beacon_sub_read_num"] = 0
            b_json["beacon_indel_percentage"] = format(100 * b_json["beacon_indel_read_num"] / b_json["beacon_aligned_read_num"], ".1f")
            b_json["beacon_sub_percentage"] = format(100 * b_json["beacon_sub_read_num"] / b_json["beacon_aligned_read_num"], ".1f")
            b_json["perfect_beacon_percent"] = format(100 * (b_json["beacon_aligned_read_num"] - b_json["beacon_indel_read_num"]) / (
                    b_json["wt_aligned_read_num"] + b_json["beacon_aligned_read_num"]), ".1f")
            b_json["beacon_fidelity"] = format(
                100 * (b_json["beacon_aligned_read_num"] - b_json["beacon_indel_read_num"]) / b_json["beacon_aligned_read_num"], ".1f")
        if ref_name == "Reference" and qw_name == "sg_cut":
            b_json["indel_read_num"] = int(stats["indels"])
            b_json["sub_read_num"] = int(stats["substitution"])
            b_json["indel_percentage"] = format(100 * b_json["indel_read_num"] / b_json["wt_aligned_read_num"], ".1f")

    df = pd.DataFrame(qw_stats)
    df.insert(0, "samplename", cs2_info["running_info"]["args"].name)
    df.to_csv(cs2_folder + "/CRISPResso_qw_stats.txt", sep="\t", header=True, index=False, na_rep=0)
    return b_json


def find_fastq_files(path, pattern):
    # Get the absolute path of the directory
    abs_path = os.path.abspath(path)
    
    # Join the pattern to the path
    file_path = os.path.join(abs_path, pattern)
    
    # Use glob to find all matches and return them
    return glob.glob(file_path, recursive=True)


def get_subdirectories(directory):
    subdirectories = []
    for entry in os.scandir(directory):
        if entry.is_dir():
            subdirectories.append(entry.path)
    return subdirectories

# i think its the assay bed file
def generate_amplicon_list(bed_file):

    genome_build = 'hg38'

    reference_path = f'/data/references/{genome_build}.fa'
    on_target_amplicon = 'CAGACGAGTCGGATCTCCCTTTGGGCCGCCTCCCCGCCTGGTACCGGCCGGCTTGTCGACGACGGCGGTCTCCGTCGTCAGGATCATCCGGTTTAAGACCAATGACTTACAAGGCAGCTGTAGATACCGTTCGT'

    with open(bed_file,'r') as f, open('amplicon_info.txt','w') as fw:
        for i,line in enumerate(f):
            curr_line = line.strip().split(',')
            # skip header
            if i == 0:
                continue
            site_name,chr,start,end,info,di_start,di_end,strand = curr_line[0],curr_line[1],curr_line[2],curr_line[3],curr_line[4],curr_line[5],curr_line[6],curr_line[7]

            amplicon_sequence_command = f'samtools faidx {reference_path} {chr}:{start}-{end}'

            amplicon_sequence = ''.join(subprocess.check_output(amplicon_sequence_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:]).upper()

            if strand == '-':
                amplicon_sequence = reverse_complement(amplicon_sequence)

            fw.write(site_name + '\t' + amplicon_sequence+'\n')
    
        fw.write('on_target' + '\t' + on_target_amplicon)
    
def generate_quant_windows(bed_file, window_size):
    quant_windows = []
    with open(bed_file,'r') as f:
        for i,line in enumerate(f):
            # skip header
            if i == 0:
                continue
            curr_line = line.strip().split(',')
            site_name,start,end,di_start,di_end,strand = curr_line[0],int(curr_line[2]),int(curr_line[3]),int(curr_line[5]),int(curr_line[6]),curr_line[7]
            seq_length = end - start
            start_relative_to_amp = di_start - start
            end_relative_to_amp = di_end - start
            if strand == '+':
                lower_bound = start_relative_to_amp - window_size + 1
                upper_bound = end_relative_to_amp + window_size + 1
            else:  # strand == '-'
                # Calculate start and end sites relative to genome differently for '-' strand.
                lower_bound = seq_length - end_relative_to_amp - window_size + 1
                upper_bound = seq_length - start_relative_to_amp + window_size + 1 
            quant_windows.append(site_name + ":dinuc:" + str(lower_bound) + "-" + str(upper_bound) + ":0")
    # on target, might be one-off for this project
    quant_windows.append('on_target:dinuc:65-72:0')
    return(quant_windows)


def run_crispresso2(sample_name,read_parent_loc,quant_windows, is_mixed):

    reads_loc = read_parent_loc + sample_name
    name = sample_name + '_results'
    
    fastq1 = find_fastq_files(reads_loc, '*R1*.fastq.gz')[0]
    fastq2 = find_fastq_files(reads_loc, '*R2*.fastq.gz')[0]
    nproc = 16

    # run pooled and mixed mode
    if is_mixed:
        genome_build = 'hg38'
        reference_path = f'/data/references/{genome_build}'
        command = "CRISPRessoPooled -r1 %s -r2 %s -x %s -f amplicon_info.txt --name %s \
            --min_frequency_alleles_around_cut_to_plot 0.05 --write_detailed_allele_table \
            --needleman_wunsch_gap_extend 0 --trim_sequences --trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:10:0:true \
            --place_report_in_output_folder --n_processes %s --bam_output" % (fastq1, fastq2, reference_path, name, nproc)
    else:
       # run pooled only
        command = "CRISPRessoPooled -r1 %s -r2 %s -f amplicon_info.txt --name %s \
            --min_frequency_alleles_around_cut_to_plot 0.05 --write_detailed_allele_table \
            --needleman_wunsch_gap_extend 0 --trim_sequences --trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:10:0:true \
            --place_report_in_output_folder --n_processes %s --bam_output" % (fastq1, fastq2, name, nproc)
   
    subprocess.call(command , shell=True)
    
    # Create an empty DataFrame to store the results
    results_df = pd.DataFrame(columns=["samplename", "amplicon", "window_name", "window_region","initial_reads","reads_after_merge","reads_align","align_perc","unmodified","modified", "total","indels", "insertion", "deletion", "substitution", "whole_window_deletion", "indel_perc"])

    # Create an empty list to store the results
    results = []

    for i, quant_window in enumerate(quant_windows):
        amplicon_name = quant_window.split(':')[0]
        path = 'CRISPRessoPooled_on_' + name + '/CRISPResso_on_' + amplicon_name

        path = path.replace('.','_')

        os.makedirs(os.path.join(path, "cs2_alignment_html"), exist_ok=True)

        allele2html_command = "python /data/tbOnT/utils/allele2html.py -f %s -r %s -b %s" % (path+'/', 'Reference', quant_window)

        subprocess.call(allele2html_command, shell=True)
        
        mapping_stats_path = path + '/CRISPResso_mapping_statistics.txt'

        if os.path.exists(mapping_stats_path):
            mapping_stats = pd.read_csv(mapping_stats_path, sep='\t')
            initial_reads = int(mapping_stats['READS IN INPUTS'].iloc[0])
            reads_after_merge = int(mapping_stats['READS AFTER PREPROCESSING'].iloc[0])
            reads_aligned = int(mapping_stats['READS ALIGNED'].iloc[0])
            align_perc = float(reads_aligned/reads_after_merge)
        else:
            # if doesn't exist, there's a mapping error. Still get num reads and merged reads from different files.

            r1_file = path + '/out.notCombined_1.fastq.gz'
            command = f'zcat {r1_file} | wc -l'
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            not_combined1_res, _ = process.communicate()
            not_combined1_reads = int(not_combined1_res) // 4

            r2_file = path + '/out.notCombined_2.fastq.gz'
            command = f'zcat {r2_file} | wc -l'
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            not_combined2_res, _ = process.communicate()
            not_combined2_reads = int(not_combined2_res) // 4

            r1r2_file = path + '/out.extendedFrags.fastq.gz'
            command = f'zcat {r1r2_file} | wc -l'
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, _ = process.communicate()
            reads_after_merge = int(output) // 4
            initial_reads = not_combined1_reads + not_combined2_reads + reads_after_merge
            reads_aligned = 0
            align_perc = 0
                
        if not os.path.exists(path):
            continue

        window_quantification(path, [quant_window])

        qw_stats_fn = path + '/CRISPResso_qw_stats.txt'

        if not os.path.exists(qw_stats_fn):
            continue

        with open(qw_stats_fn, 'r') as f:
            line = f.readlines()[1]

        if len(line.strip().split('\t')) == 11:
            samplename, amplicon, window_name, window_region, modified, indels, insertion, deletion, substitution, whole_window_deletion, unmodified = line.strip().split('\t')
        else:
            continue
        total = int(modified) + int(unmodified)
        indel_perc = np.round(100 * int(indels) / (total), 6)
            
        # Append the results to the list
        results.append({
            "samplename": samplename,
            "amplicon": amplicon,
            "window_name": window_name,
            "window_region": window_region,
            "initial_reads": initial_reads,
            "reads_after_merge": reads_after_merge,
            "reads_aligned": reads_aligned,
            "align_perc": align_perc,
            "unmodified": unmodified,
            "modified": modified,
            "total": total,
            "indels": indels,
            "insertion": insertion,
            "deletion": deletion,
            "substitution": substitution,
            "whole_window_deletion": whole_window_deletion,
            "indel_perc": indel_perc
        })
    
    # Create a DataFrame from the results list
    results_df = pd.DataFrame(results)

    # Modify column names after the fourth column

    for column in results_df.columns[4:]:
        results_df.rename(columns={column: sample_name + '_' + column}, inplace=True)
        
    return(results_df)

def merge_dfs(df_list):
    # Join DataFrames based on the values in the first four columns
    merged_df = df_list[0]  # Initialize with the first DataFrame

    for df in df_list[1:]:
        merged_df = pd.merge(merged_df, df, on=["samplename", "amplicon", "window_name", "window_region"], how="outer")

    # Print the merged DataFrame
    return(merged_df)

if __name__ == '__main__':

    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description='Process bed_file and sample_name arguments.')

    # Add arguments
    parser.add_argument('--project-name', dest='project_name',
                        help='Name of project')
    parser.add_argument('--bed-file', dest='bed_file',
                        help='Path to the bed file')
    parser.add_argument('--window-size', dest='window_size',
                        help='Window size around region of interest')
    parser.add_argument('--reads', dest='reads',
                        help='Parent directory where all reads are stored')
    parser.add_argument('--mixed', action='store_true', help='Enable mixed mode (aligns to genome and amplicon)')
    # parser.add_argument('--sample-names', dest='sample_names', nargs='+', 
    #                     help='Sample names')


    # Parse the command-line arguments
    args = parser.parse_args()

    # Access the values using the argument names
    bed_file = args.bed_file
    #sample_names = args.sample_names
    window_size = int(args.window_size)
    project_name = args.project_name
    reads_dir = args.reads
    is_mixed = args.mixed

    sample_names = []
    # Example usage
    subdirectories = get_subdirectories(reads_dir)
    for subdirectory in subdirectories:
        sample_names.append(subdirectory.split('/')[-1])
    
    generate_amplicon_list(bed_file)
    dfs = []
    quant_windows = generate_quant_windows(bed_file, window_size)
    for sample_num, sample_name in enumerate(sample_names):
        print('Running CRISPResso2 on ' + sample_name + '...')
        df = run_crispresso2(sample_name,reads_dir,quant_windows,is_mixed)
        dfs.append(df)
    merged_df = merge_dfs(dfs)
    merged_df.to_csv(project_name+"_results.csv", index=False)

