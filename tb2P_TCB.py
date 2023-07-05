#!/usr/bin/env python
"""
Tome 2P-seq analysis
Input from Benchling
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
from io import StringIO
import re
import quilt3
from benchling_sdk.benchling import Benchling
from benchling_sdk.auth.api_key_auth import ApiKeyAuth

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio  # free up some memory
        sys.stdout = self._stdout

def download_sequencing_data(genomics_ampseq_project_queue,reads_parent_dir):
    # get list of projects
    raw_project_list = subprocess.check_output('bs list projects -f csv',shell=True).decode(sys.stdout.encoding)
    raw_project_list = io.StringIO(raw_project_list)

    for line in raw_project_list:
        if genomics_ampseq_project_queue in line:
            project_id = line.split(',')[1]
    raw_project_list.close()

    download_folder = reads_parent_dir + genomics_ampseq_project_queue
    download_fastq_command = 'bs download project -i ' + str(project_id) + ' -o ' + download_folder + ' --extension=fastq.gz'
    subprocess.check_output(download_fastq_command, shell=True)
    subdirectories = [subdirectory for subdirectory in os.listdir(download_folder) if os.path.isdir(os.path.join(download_folder, subdirectory))]
    return(subdirectories)

def find_fastq_files(path, pattern):
    # Get the absolute path of the directory
    abs_path = os.path.abspath(path)
    
    # Join the pattern to the path
    file_path = os.path.join(abs_path, pattern)
    
    # Use glob to find all matches and return them
    return glob.glob(file_path, recursive=True)

def find_directory(directory, substring):
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            if substring in dir:
                return os.path.join(root, dir)

    return None  # Return None if no directory was found.

def download_plasmid(payload_id):      
    query = '''
    SELECT 
        url$ 
    FROM
        plasmid
    WHERE
        file_registry_id$ = %s
    '''
    cur.execute(query, [payload_id])
    result = cur.fetchone()    
    if result:
        url = result[0]
        pattern = r"(seq_[A-Za-z0-9]+)"
        match = re.search(pattern, url)
        if match:
            seq_id = match.group(1)
            print(seq_id)
    benchling = Benchling(url=api_url, auth_method=ApiKeyAuth(api_key))
    
    sequence_generator = benchling.dna_sequences.list()

    found_sequence = False
    for page in sequence_generator:
        for sequence in page:
            if sequence.id == seq_id:
                bases = sequence.bases.upper()
                found_sequence = True
                break
        if found_sequence:
            break
    output_fn = f'/data/references/{payload_id}.fa'
    print(output_fn)
    with open(output_fn,'w') as fw:
        fw.write(f'>{payload_id}\n{bases}')

def generate_amplicon(genomics_ampseq_project_queue):

    amplicon_dict = {}

    # Single SQL query to get all necessary data
    query = '''
        SELECT 
            basespace_sample_name, 
            atg_atg.file_registry_id$,
            payload.file_registry_id$,
            p1.chromosome,
            p1.start,
            p1.end,
            p1.adapter_coordinates,
            p1.genome_build,
            p1.strand,
            dna_oligo3.bases,
            p2.start,
            p2.end,
            p2.adapter_coordinates,
            p2.strand,
            dna_oligo4.bases,
            dna_oligo1.bases,
            spacer1.cut_position,
            spacer1.genome_build,
            spacer1.strand,
            dna_oligo2.bases,
            spacer2.cut_position,
            spacer2.genome_build,
            spacer2.strand,
            beacon1.bases as beacon1_bases,
            beacon2.bases as beacon2_bases,
            attb.left_half,
            attb.right_half,
            attp.left_half,
            attp.right_half,
            attb_site.bases as attb_bases,
            attp_site.bases as attp_bases,
            target_gene.direction_of_transcription
        FROM 
            junction_sequencing_sample_metasheet$raw
            JOIN atg_atg AS atg_atg ON guide_id = atg_atg.id 
            JOIN plasmid AS payload ON payload_id = payload.id
            JOIN primer_pair as pp on primer_pair_id = pp.id
            JOIN primer as p1 on pp.forward_primer = p1.id
            JOIN primer as p2 on pp.reverse_primer = p2.id
            JOIN modified_rna AS mod_rna1 ON mod_rna1.id = atg_atg.atg1
            JOIN modified_rna AS mod_rna2 ON mod_rna2.id = atg_atg.atg2
            JOIN atgRNA AS atg_rna1 ON atg_rna1.id = mod_rna1.rna
            JOIN atgRNA AS atg_rna2 ON atg_rna2.id = mod_rna2.rna
            JOIN spacer as spacer1 ON spacer1.id = atg_rna1.spacer
            JOIN spacer as spacer2 ON spacer2.id = atg_rna2.spacer
            JOIN dna_oligo AS dna_oligo1 ON dna_oligo1.id = spacer1.id
            JOIN dna_oligo AS dna_oligo2 ON dna_oligo2.id = spacer2.id 
            JOIN dna_oligo AS dna_oligo3 ON dna_oligo3.id = p1.id
            JOIN dna_oligo AS dna_oligo4 ON dna_oligo4.id = p2.id
            JOIN dna_sequence as beacon1 on beacon1.id = atg_rna1.beacon 
            JOIN dna_sequence as beacon2 on beacon2.id = atg_rna2.beacon 
            JOIN attachment_sequence as attb on attb.id = atg_atg.expected_beacon
            JOIN attachment_sequence as attp on attp.id = payload.attachment_sequence
            JOIN dna_sequence as attb_site on attb_site.id = atg_atg.expected_beacon
            JOIN dna_sequence as attp_site on attp_site.id = payload.attachment_sequence
            JOIN target_gene as target_gene ON target_gene.id = spacer1.target_gene
        WHERE genomics_sequencing_project_queue = %s
        '''
    cur.execute(query, [genomics_ampseq_project_queue])   
    #result = cur.fetchone()    
    for result in cur.fetchall():
        if result:
            bs_sample_name, aaanpnsg_id, payload_id, gp_chr, gp_start, gp_end, gp_adapter_coords, gp_build, gp_strand, gp_bases, cp_start, cp_end, cp_adapter_coords, cp_strand, cp_bases, spacer1_bases, spacer1_cut_pos, spacer1_genome_build, spacer1_strand, spacer2_bases, spacer2_cut_pos, spacer2_genome_build, spacer2_strand, beacon1_seq, beacon2_seq, attb_seq_left_coords, attb_seq_right_coords, attp_seq_left_coords, attp_seq_right_coords,attb_seq, attp_seq,gene_strand = result
            print(bs_sample_name)
            attb_reg_coords = [int(coord.replace('[','').replace(']','')) for coord in attb_seq_left_coords.split(',')]
            attb_prime_coords = [int(coord.replace('[','').replace(']','')) for coord in attb_seq_right_coords.split(',')]
            attp_reg_coords = [int(coord.replace('[','').replace(']','')) for coord in attp_seq_left_coords.split(',')]
            attp_prime_coords = [int(coord.replace('[','').replace(']','')) for coord in attp_seq_right_coords.split(',')]

            ATTB_REG = attb_seq[attb_reg_coords[0]-1:attb_reg_coords[1]+2].upper()
            ATTB_PRIME = attb_seq[attb_prime_coords[0]-1:attb_prime_coords[1]-1].upper()

            ATTP_REG = attp_seq[attp_reg_coords[0]-1:attp_reg_coords[1]+2].upper()
            ATTP_PRIME = attp_seq[attp_prime_coords[0]-1:attp_prime_coords[1]-1].upper()

            cp_adapter_coords = [int(coord.replace('[','').replace(']','')) for coord in cp_adapter_coords.split(',')]

            gp_adapter_coords = [int(coord.replace('[','').replace(']','')) for coord in gp_adapter_coords.split(',')]
            
            if spacer1_cut_pos > spacer2_cut_pos:
                spacer1_cut_pos, spacer2_cut_pos = spacer2_cut_pos, spacer1_cut_pos
                spacer1_strand, spacer2_strand = spacer2_strand, spacer1_strand
        
            reference_path = f'/data/references/{gp_build}.fa'

            if gene_strand == 'antisense' or gene_strand == '-':
                # AttL
                if gp_start > spacer2_cut_pos:
                    beginning_sequence_command = f'samtools faidx {reference_path} {gp_chr}:{spacer2_cut_pos}-{gp_end}'
                # AttR
                else:
                    beginning_sequence_command = f'samtools faidx {reference_path} {gp_chr}:{gp_start}-{spacer1_cut_pos-1}'
            else:
                # AttR
                if gp_start > spacer2_cut_pos:
                    beginning_sequence_command = f'samtools faidx {reference_path} {gp_chr}:{spacer2_cut_pos}-{gp_end}'
                # AttL
                else:
                    beginning_sequence_command = f'samtools faidx {reference_path} {gp_chr}:{gp_start}-{spacer1_cut_pos-1}'

            beginning_sequence = ''.join(subprocess.check_output(beginning_sequence_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:])

            cargo_path = f'/data/references/{payload_id}.fa'

            if not os.path.exists(cargo_path):
                download_plasmid(payload_id)

            # store cargo sequence so we can search for subsequence locations
            with open(cargo_path, 'r') as f:
                cargo_sequence = ''.join(line.strip() for line in f if not line.startswith('>')).upper()

            cargo_primer_sequence_without_illumina_adapter = cp_bases[cp_adapter_coords[-1]:]
            cp_end = cargo_sequence.find(reverse_complement(cargo_primer_sequence_without_illumina_adapter)) + len(cargo_primer_sequence_without_illumina_adapter)

            # end of P' loc
            attp_prime_end  = cargo_sequence.find(ATTP_PRIME) + len(ATTP_PRIME)
            attp_reg_start = cargo_sequence.find(ATTP_REG)
            
            # proxy for direction of transcription? 
            if gene_strand == 'antisense' or gene_strand == '-':
                beginning_sequence = reverse_complement(beginning_sequence)

            # AttL
            if gene_strand == 'sense' or gene_strand == '+':
                if gp_start < spacer1_cut_pos:
                    complete_amplicon = beginning_sequence + ATTB_REG + ATTP_PRIME + cargo_sequence[attp_prime_end:cp_end]
                    sequence_to_highlight = ATTB_REG + ATTP_PRIME
                else:
                    complete_amplicon = cargo_sequence[cp_start:attp_reg_start] + ATTP_REG + ATTB_PRIME + beginning_sequence
                    sequence_to_highlight = ATTP_REG + ATTB_PRIME
            elif gene_strand == 'antisense' or gene_strand == '-':
                if gp_start > spacer1_cut_pos:
                    complete_amplicon = beginning_sequence + ATTB_REG + ATTP_PRIME + cargo_sequence[attp_prime_end:cp_end]
                    sequence_to_highlight = ATTB_REG + ATTP_PRIME
                else:
                    complete_amplicon = cargo_sequence[cp_start:attp_reg_start] + ATTP_REG + ATTB_PRIME + beginning_sequence
                    sequence_to_highlight = ATTP_REG + ATTB_PRIME

            # print(f'beginning sequence: {beginning_sequence}')
            # print(f'cargo sequence: {cargo_sequence}')
            # print(f'cargo primer part: {cargo_sequence[attp_reg_start:cp_end]}')
            # print(f'B reg: {ATTB_REG}')
            # print(f'B prime: {ATTB_PRIME}')
            # print(f'P reg: {ATTP_REG}')
            # print(f'P prime: {ATTP_PRIME}')
            print(f'Amplicon: {complete_amplicon}')
            print('Length of amplicon: ' + str(len(complete_amplicon)))
            amplicon_dict[bs_sample_name] = (complete_amplicon,gene_strand,sequence_to_highlight)
    return(amplicon_dict)

def run_crispresso2(complete_amplicon, fastq1, fastq2, name, output_folder, n_processes, gene_strand,genomics_ampseq_project_queue,quant_window_string):
    
    if gene_strand == 'sense' or gene_strand == '+':
        subprocess.call(
                "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name Cargo --min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s --write_detailed_allele_table --place_report_in_output_folder --trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:10:0:true --bam_output --n_processes %s" % (
                    fastq1, fastq2, complete_amplicon, name, output_folder, n_processes), shell=True)
        mapping_stats_path = './' + output_folder + '/CRISPResso_on_' + name + '/CRISPResso_mapping_statistics.txt'
        print(mapping_stats_path)
        if os.path.exists(mapping_stats_path):
            mapping_stats = pd.read_csv(mapping_stats_path, sep='\t')
            print(mapping_stats)
            initial_reads = int(mapping_stats['READS IN INPUTS'])
            reads_after_merge = int(mapping_stats['READS AFTER PREPROCESSING'])
            reads_aligned_to_cargo = int(mapping_stats['READS ALIGNED'])
            align_perc = float(reads_aligned_to_cargo/reads_after_merge)
            print(align_perc)
        else:
            # if doesn't exist, there's a mapping error. Still get num reads and merged reads from different files.

            r1_file = '/data/tbOnT/' + output_folder + '/CRISPResso_on_' + name + '/out.notCombined_1.fastq.gz'
            command = f'zcat {r1_file} | wc -l'
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            not_combined1_res, _ = process.communicate()
            not_combined1_reads = int(not_combined1_res) // 4

            r2_file = '/data/tbOnT/' + output_folder + '/CRISPResso_on_' + name + '/out.notCombined_2.fastq.gz'
            command = f'zcat {r2_file} | wc -l'
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            not_combined2_res, _ = process.communicate()
            not_combined2_reads = int(not_combined2_res) // 4

            r1r2_file = '/data/tbOnT/' + output_folder + '/CRISPResso_on_' + name + '/out.extendedFrags.fastq.gz'
            command = f'zcat {r1r2_file} | wc -l'
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, _ = process.communicate()
            reads_after_merge = int(output) // 4
            initial_reads = not_combined1_reads + not_combined2_reads + reads_after_merge
            reads_aligned_to_cargo = 0
            align_perc = 0
    else:
        subprocess.call(
                "CRISPResso --fastq_r2 %s --fastq_r1 %s --amplicon_seq %s --amplicon_name Cargo --min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s --write_detailed_allele_table --place_report_in_output_folder --trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:10:0:true --bam_output --n_processes %s" % (
                    fastq1, fastq2, complete_amplicon, name, output_folder, n_processes), shell=True)
        mapping_stats_path = './' + output_folder + '/CRISPResso_on_' + name + '/CRISPResso_mapping_statistics.txt'
        print(mapping_stats_path)
        if os.path.exists(mapping_stats_path):
            mapping_stats = pd.read_csv(mapping_stats_path, sep='\t')
            print(mapping_stats)
            initial_reads = int(mapping_stats['READS IN INPUTS'])
            reads_after_merge = int(mapping_stats['READS AFTER PREPROCESSING'])
            reads_aligned_to_cargo = int(mapping_stats['READS ALIGNED'])
            align_perc = float(reads_aligned_to_cargo/reads_after_merge)
        else:
            # if doesn't exist, there's a mapping error. Still get num reads and merged reads from different files.

            r1_file = '/data/tbOnT/' + output_folder + '/CRISPResso_on_' + name + '/out.notCombined_1.fastq.gz'
            command = f'zcat {r1_file} | wc -l'
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            not_combined1_res, _ = process.communicate()
            not_combined1_reads = int(not_combined1_res) // 4

            r2_file = '/data/tbOnT/' + output_folder + '/CRISPResso_on_' + name + '/out.notCombined_2.fastq.gz'
            command = f'zcat {r2_file} | wc -l'
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            not_combined2_res, _ = process.communicate()
            not_combined2_reads = int(not_combined2_res) // 4

            r1r2_file = '/data/tbOnT/' + output_folder + '/CRISPResso_on_' + name + '/out.extendedFrags.fastq.gz'
            command = f'zcat {r1r2_file} | wc -l'
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, _ = process.communicate()
            reads_after_merge = int(output) // 4
            initial_reads = not_combined1_reads + not_combined2_reads + reads_after_merge
            reads_aligned_to_cargo = 0
            align_perc = 0

    curr_folder_path = name + '/' + "CRISPResso_on_" + name
    os.makedirs(curr_folder_path + '/cs2_alignment_html', exist_ok=True)
    allele2html_command = "python /data/tbOnT/utils/allele2html.py -f %s -r %s -b %s" % (curr_folder_path+'/', 'Cargo', quant_window_string)
    #allele2html_command = "python /data/tbOnT/utils/allele2html.py -f %s -r %s" % (curr_folder_path+'/', 'Cargo')
    subprocess.call(allele2html_command, shell=True)
    return(initial_reads,reads_after_merge,reads_aligned_to_cargo,align_perc)
    
if __name__ == "__main__":
    # # Create the parser
    # parser = argparse.ArgumentParser(description='This script accepts a genomics project queue as an argument.')
    # # Add an argument
    # parser.add_argument('--project_name', type=str, required=True, help='Project name as indicated in ampseq_genomics_project_queue column.')
    # # Parse the arguments
    # args = parser.parse_args()
    # genomics_ampseq_project_queue = args.project_name
    # reads_parent_dir = '/data/reads/'
    # ncpu = 8
    # new_data_subdirs = download_sequencing_data(genomics_ampseq_project_queue,reads_parent_dir)
    # print('Generating expected amplicon for ' + genomics_ampseq_project_queue)
    # amplicon_dict = generate_amplicon(genomics_ampseq_project_queue)
    # sample_names = []
    # num_initial_reads = []
    # num_reads_after_merge = []
    # num_reads_aligned_to_cargo = []
    # align_percs = []
    # print(amplicon_dict)
    # for sample_name in amplicon_dict:
    #     curr_amplicon = amplicon_dict[sample_name][0]
    #     gene_strand = amplicon_dict[sample_name][1]
    #     sequence_to_highlight = amplicon_dict[sample_name][2]
    #     print(sequence_to_highlight)
    #     quant_window_lower = curr_amplicon.find(sequence_to_highlight) + 1
    #     quant_window_upper = quant_window_lower + len(sequence_to_highlight) - 1
    #     quant_window_string = "on_target:junction:" + str(quant_window_lower) + "-" + str(quant_window_upper) + ":0"
    #     reads_project_dir = f'{reads_parent_dir}{genomics_ampseq_project_queue}'
    #     reads_loc = find_directory(reads_project_dir, sample_name)
    #     fastq1 = find_fastq_files(reads_loc, '*_R1_*.fastq.gz')[0]
    #     fastq2 = find_fastq_files(reads_loc, '*_R2_*.fastq.gz')[0]
    #     print('Running cs2 for ' + sample_name + '...')
    #     initial_reads, reads_after_merge, reads_aligned_to_cargo, align_perc = run_crispresso2(curr_amplicon, fastq1, fastq2,sample_name,sample_name,ncpu,gene_strand,genomics_ampseq_project_queue, quant_window_string)
    #     # Append the results to the lists
    #     sample_names.append(sample_name)
    #     num_initial_reads.append(initial_reads)
    #     num_reads_after_merge.append(reads_after_merge)
    #     num_reads_aligned_to_cargo.append(reads_aligned_to_cargo)
    #     align_percs.append(round(100*align_perc,2))
    
    # assert len(num_initial_reads) == len(num_reads_after_merge) == len(num_reads_aligned_to_cargo) == len(sample_names) == len(align_percs)

    # # Specify the file name
    # output_fn = f'{genomics_ampseq_project_queue}_data.csv'

    # # Open the file in write mode
    # with open(output_fn, 'w') as f:
    #     # Write the header row
    #     header = 'Sample name,Total reads,Reads after merging R1/R2,Reads aligning to Cargo,Align %\n'
    #     f.write(header)

    #     # Write each row
    #     for i in range(len(num_initial_reads)):
    #         row = f'{sample_names[i]},{num_initial_reads[i]},{num_reads_after_merge[i]},{num_reads_aligned_to_cargo[i]},{align_percs[i]}\n'
    #         f.write(row)

    # print(f'Data written to {output_fn} successfully.')

    # check if the package existed
    ngs_id = "BTB201"
    pipeline_run_id = "MJK_2PrimerJunction_BTB201a_20230609"

    if "AmpSeq/" + ngs_id in list(quilt3.list_packages("s3://tb-ngs-quilt/")):
        quilt3.Package.install("AmpSeq/" + ngs_id, "s3://tb-ngs-quilt/")
        p = quilt3.Package.browse("AmpSeq/" + ngs_id)
    else:
        p = quilt3.Package()

    # adding data
    # input package
    # p.set_dir("fastq/" + pipeline_run_id, pipeline_run_id)

    # output package
    p.set(pipeline_run_id + "/stats.csv", "MJK_2PrimerJunction_BTB201a_20230609_data.csv")
    p.set_dir(pipeline_run_id + "/cs2_alignment_html",  "cs2_alignment_html/")

    # Pushing a package to a remote registry
    with Capturing() as output:
        p.push("AmpSeq/" + ngs_id, "s3://tb-ngs-quilt/", force=True)
    base_url = output[1].split()[-1]
    full_url = f"{base_url}/tree/{p.top_hash}"
    print(full_url)