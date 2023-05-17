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

def generate_amplicon(genomics_ampseq_project_queue):   

    # get primer, AA, payload for further SQL queries
    cur.execute('SELECT sample_name, primer1.file_registry_id$, primer2.file_registry_id$, atg_atg.file_registry_id$, payload.file_registry_id$ FROM _2pseq_sample_metasheet$raw JOIN atg_atg AS atg_atg ON aaanpnsg_id = atg_atg.id JOIN primer AS primer1 ON genomic_primer = primer1.id JOIN primer AS primer2 ON cargo_primer = primer2.id JOIN plasmid AS payload ON payload_id = payload.id WHERE genomics_ampseq_project_queue = %s',[genomics_ampseq_project_queue])
    sample_name, genomic_primer, cargo_primer, aaanpnsg_id, payload_id = cur.fetchone()

    # retrieve info for genomic primer
    cur.execute("select chromosome, start, p1.end, adapter_coordinates, genome_build, strand from primer as p1 where file_registry_id$ = %s", [genomic_primer])
    gp_chr, gp_start, gp_end, gp_adapter_coords, gp_build, gp_strand = cur.fetchone()

    # retrieve info for cargo primer
    cur.execute("select start, p1.end, p1.adapter_coordinates, p1.strand from primer as p1 where file_registry_id$ = %s", [cargo_primer])
    cp_start, cp_end, cp_adapter_coords, cp_strand = cur.fetchone()

    cp_adapter_lower_bound = int(cp_adapter_coords.split(',')[0].split('[')[-1])
    cp_adapter_upper_bound = int(cp_adapter_coords.split(',')[1].split(']')[0])
    cp_adapter_coords = [cp_adapter_lower_bound, cp_adapter_upper_bound]

    gp_adapter_lower_bound = int(gp_adapter_coords.split(',')[0].split('[')[-1])
    gp_adapter_upper_bound = int(gp_adapter_coords.split(',')[1].split(']')[0])
    gp_adapter_coords = [gp_adapter_lower_bound, gp_adapter_upper_bound]

    # get cut position, genome build, strand for both atg1 spacer and atg2 spacer
    # change "sequence" to "rna" once on main site
    cur.execute("SELECT dna_oligo1.bases, spacer1.cut_position, spacer1.genome_build, spacer1.strand, dna_oligo2.bases, spacer2.cut_position, spacer2.genome_build,spacer2.strand FROM atg_atg JOIN modified_rna AS mod_rna1 ON mod_rna1.id = atg_atg.atg1 JOIN modified_rna AS mod_rna2 ON mod_rna2.id = atg_atg.atg2 JOIN atgRNA AS atg_rna1 ON atg_rna1.id = mod_rna1.rna JOIN atgRNA AS atg_rna2 ON atg_rna2.id = mod_rna2.rna JOIN spacer as spacer1 ON spacer1.id = atg_rna1.spacer JOIN spacer as spacer2 ON spacer2.id = atg_rna2.spacer JOIN dna_oligo AS dna_oligo1 ON dna_oligo1.id = spacer1.id JOIN dna_oligo AS dna_oligo2 ON dna_oligo2.id = spacer2.id WHERE atg_atg.file_registry_id$ = %s",[aaanpnsg_id])
    spacer1_bases, spacer1_cut_pos, spacer1_genome_build, spacer1_strand, spacer2_bases, spacer2_cut_pos, spacer2_genome_build, spacer2_strand = cur.fetchone()

    # get genomic sequence from beginning of genomic primer to spacer 1 cut site 
    reference_path = '/data/references/' + gp_build + '.fa'

    # get atg1 beacon and atg2 beacon sequences
    # change "sequence" to "rna" once on main site
    cur.execute("SELECT beacon1.bases as beacon1_bases, beacon2.bases as beacon2_bases from atg_atg JOIN modified_rna AS mod_rna1 ON atg1 = mod_rna1.id JOIN modified_rna AS mod_rna2 ON atg2 = mod_rna2.id JOIN atgRNA AS atg_rna1 ON mod_rna1.rna = atg_rna1.id JOIN atgRNA AS atg_rna2 ON mod_rna2.rna = atg_rna2.id JOIN dna_sequence as beacon1 on beacon1.id = atg_rna1.beacon JOIN dna_sequence as beacon2 on beacon2.id = atg_rna2.beacon WHERE atg_atg.file_registry_id$ = %s",[aaanpnsg_id])
    beacon1_seq, beacon2_seq = cur.fetchone()

    
    beacon_overlap_sequence = get_beacon_seq(beacon1_seq, spacer1_strand, beacon2_seq, spacer2_strand)

    if spacer1_cut_pos > spacer2_cut_pos:
        beacon_overlap_sequence = reverse_complement(beacon_overlap_sequence)
        beginning_sequence_command = 'samtools faidx ' +  reference_path + ' ' + gp_chr + ':' + str(gp_start) + '-' + str(spacer2_cut_pos-1)
        # run samtools faidx, get output, skip header line, concatenate rest of lines to get full sequence
        beginning_sequence = ''.join(subprocess.check_output(beginning_sequence_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:])
    else:
        beginning_sequence_command = 'samtools faidx ' +  reference_path + ' ' + gp_chr + ':' + str(gp_start) + '-' + str(spacer1_cut_pos)

        # run samtools faidx, get output, skip header line, concatenate rest of lines to get full sequence
        beginning_sequence = ''.join(subprocess.check_output(beginning_sequence_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:])

    # first part of PASTE -- insert beacon sequence 
    full_sequence = beginning_sequence + beacon_overlap_sequence
    
    attp_reg = 'GTGGTTTGTCTGGTCAACCACCGCGAC'
    attp_prime = 'CTCAGTGGTGTACGGTACAAACCCA'

    attb_reg = 'GGCTTGTCGACGACGGCGAC'
    attb_prime = 'CTCCGTCGTCAGGATCAT'
    

    # path to plasmid containing cargo
    cargo_path = '/data/references/' + payload_id + '.fa'

    # store cargo sequence so we can search for subsequence locations
    cargo_sequence = ''
    with open(cargo_path, 'r') as f:
        for line in f:
            if not line.startswith('>'): # ignore header line
                cargo_sequence += line.strip()
    cargo_sequence = cargo_sequence.upper()
    
    if cp_start is None or cp_end is None:
        cur.execute('SELECT DISTINCT(bases) FROM _2pseq_sample_metasheet$raw JOIN dna_oligo as dna_oligo1 ON cargo_primer = dna_oligo1.id WHERE genomics_ampseq_project_queue = %s',[genomics_ampseq_project_queue])
        cargo_primer_sequence = cur.fetchone()[0]
        #print(cargo_primer_sequence)
        adapter_start = cp_adapter_coords[0]
        adapter_end = cp_adapter_coords[1]
        cargo_primer_sequence_without_illumina_adapter = cargo_primer_sequence[adapter_end:]
        cp_end = cargo_sequence.find(reverse_complement(cargo_primer_sequence_without_illumina_adapter)) + len(cargo_primer_sequence_without_illumina_adapter)
    
    #testing purposes only since initial 2p experiment didn't work -- this is a fake primer that exists near AttR boundary 
    cp_start = 157
    cp_end = 195

    # end of P' loc
    attp_prime_end  = cargo_sequence.find(attp_prime) + len(attp_prime)

    attp_reg_start = cargo_sequence.find(attp_reg)

    if spacer1_cut_pos > spacer2_cut_pos:
        complete_amplicon = cargo_sequence[cp_start:attp_reg_start] + (attp_reg) + (attb_prime)  + reverse_complement(beginning_sequence)
    else:
        complete_amplicon = beginning_sequence + attb_reg + attp_prime + cargo_sequence[attp_prime_end:cp_end]
    
    
    print('beginning sequence: ' + beginning_sequence)
    print('cargo sequence: ' + cargo_sequence)
    print('cargo primer part: ' + cargo_sequence[cp_start:attp_reg_start])
    print('B prime: ' + attb_prime)
    print('P reg: ' + attp_reg)
    print('Amplicon: ' + complete_amplicon)
    return(complete_amplicon)

def run_crispresso2(complete_amplicon, fastq1, fastq2, name, output_folder, n_processes, genomics_ampseq_project_queue):
    subprocess.call(
            "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name Cargo --min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s --write_detailed_allele_table --place_report_in_output_folder --n_processes %s" % (
                fastq1, fastq2, complete_amplicon, name, output, ncpu), shell=True)
    
    mapping_stats = pd.read_csv('~/projects/tbOnT/'+output+'/CRISPResso_on_'+name+'/CRISPResso_mapping_statistics.txt',sep='\t')

    print(mapping_stats)
    initial_reads = int(mapping_stats['READS IN INPUTS'])
    reads_after_merge = int(mapping_stats['READS AFTER PREPROCESSING'])
    reads_aligned_to_cargo = int(mapping_stats['READS ALIGNED'])
    align_perc = float(reads_aligned_to_cargo/reads_after_merge)
    
    print('Total reads' + '\t' + 'Reads after merging R1/R2' + '\t' + 'Reads aligning to Cargo' + '\t' + 'Align %')
    print(str(initial_reads) + '\t' + str(reads_after_merge) + '\t' + str(reads_aligned_to_cargo) + '\t' + str(np.round(100*align_perc,2)) + '%')

if __name__ == "__main__":
    genomics_ampseq_project_queue = 'LH_2Primer_TB000158c_20230428'
    reads_parent_dir = '/data/reads/'
    ncpu = 8
    new_data_subdirs = download_sequencing_data(genomics_ampseq_project_queue,reads_parent_dir)
    print('Generating expected amplicon for ' + genomics_ampseq_project_queue)
    amplicon = generate_amplicon(genomics_ampseq_project_queue)
    for subdir in new_data_subdirs:
        sample_name = subdir.split('.')[0].split('_')[0]
        fastq1 = [file for file in os.listdir(directory) if os.path.isfile(os.path.join(directory, file)) and "R1" in file]
        fastq2 = [file for file in os.listdir(directory) if os.path.isfile(os.path.join(directory, file)) and "R2" in file]
        print('Running cs2 for ' + sample_name + '...')
        #run_crispresso2(amplicon, fastq1, fastq2,sample_name,sample_name,ncpu,genomics_ampseq_project_queue)

