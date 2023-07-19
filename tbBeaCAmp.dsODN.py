#!/usr/bin/env python
"""
On-Target analysis for Amp-seq dsODN
Information from Benchling
Input SP ID
Input PRIP ID
"""

import argparse
import glob
import datetime
from utils.base import *
from utils.common_functions import *


def main():
    # Parse command line options
    parser = argparse.ArgumentParser(description='Process meta data', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", help='Path to fastq')
    parser.add_argument("-g", help='SP id')
    parser.add_argument("-a", help='prime pair id')
    parser.add_argument("-s", help='specifiy a sample in the fastq path', default=None)
    parser.add_argument("-p", help='Number of CPUs to use', default=4)
    parser.add_argument("-o", help='Output folder', default="./")
    parser.add_argument("-cs2", help='CRISPRESSO2 parameters', default="")
    parser.add_argument("-dsodn", help='dsODN seq', default="")

    args = parser.parse_args()
    fastq = args.i
    name = args.s
    aaan_id = args.g
    pp_id = args.a
    ncpu = int(args.p)
    cs2 = args.cs2
    output = args.o
    dsodn = args.dsodn

    # create output folder
    os.makedirs(os.path.join(output, "cs2_alignment_html"), exist_ok=True)

    # amplicon information file
    amplicon_fh = open(os.path.join(output, name + ".amplicon.txt"), 'w')

    # Get primer information
    cur.execute("select p1.chromosome, p1.start, p2.end, p1.genome_build, target_gene.direction_of_transcription from primer_pair "
                "join primer as p1 on p1.id = primer_pair.forward_primer "
                "join primer as p2 on p2.id = primer_pair.reverse_primer "
                "join target_gene on target_gene.id = p1.gene_or_target_name "
                "where primer_pair.file_registry_id$ = %s", [pp_id])
    if cur.rowcount == 0:
        print("Can't find this PRIP in Benchling!")
        return
    else:
        chr, wt_start, wt_end, genome_build, target_strand = cur.fetchone()

    cs2_stats = {"aaanid": aaan_id, "ppid": pp_id, "samplename": name}

    # skip if no sample name or no fastq
    if not name or len(glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R1_*")) == 0:
        print("Sample is not in the fastq path")

    # reference genome
    genome_build = re.sub(".*/", "", genome_build)
    genome_fa = "/home/ubuntu/annotation/2bit/" + genome_build + ".2bit"

    # get r1 and r2 fastq
    if target_strand == "antisense" or target_strand == "-":
        r1 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R2_*")[0]
        r2 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R1_*")[0]
    else:
        r1 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R1_*")[0]
        r2 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R2_*")[0]

    # sample job log
    job_fh = open(os.path.join(output, name + ".job.log"), 'wb')

    # WT amplicon
    wt_amplicon = get_seq(genome_fa, chr, wt_start, wt_end, target_strand)
    if len(wt_amplicon) > 298:
        cs2 = "--force_merge_pairs "
    elif len(wt_amplicon) >= 293 and len(wt_amplicon) <= 298:
        cs2 = "--stringent_flash_merging "
    amplicon_fh.write(name + "\tWT\t" + wt_amplicon + "\n")

    sp1_info = {}

    # Beacon amplicon
    # control sample w/o AA/SG/AN/PN id
    if not dsodn:
        subprocess.call(
            "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT --name %s --output_folder %s --n_processes %s "
            "--min_frequency_alleles_around_cut_to_plot 0.05 --write_detailed_allele_table --needleman_wunsch_gap_extend 0 "
            "--trim_sequences --trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:10:0:true "
            "--bam_output --suppress_report --place_report_in_output_folder %s " % (r1, r2, wt_amplicon, name, output, ncpu, cs2),
            stderr=job_fh, stdout=job_fh, shell=True)

    # dsODN
    else:
        # get spacer sequences, beacon sequences, ngRNA sequences
        cur.execute("select sp.bases from spacer "
                    "join dna_oligo as sp on sp.id = spacer.id "
                    "where spacer.file_registry_id$ = %s", [aaan_id])
        sp_seq = cur.fetchone()[0]
        print(sp_seq)
        sp1_info = get_cut_site(wt_amplicon, sp_seq)

        beacon = dsodn

        # beacon seq
        beacon_amplicon = wt_amplicon[0:sp1_info["cut"]] + beacon + wt_amplicon[sp1_info["cut"]:]
        if len(beacon_amplicon) > 298:
            cs2 = "--force_merge_pairs "
        elif len(beacon_amplicon) >= 293 and len(beacon_amplicon) <= 298:
            cs2 = "--stringent_flash_merging "

        amplicon_fh.write(name + "\tBeacon\t" + beacon_amplicon + "\n")

        # define quantification window
        # WT amplicon, spacer cutting 2bp
        wt_qw1 = "WT:spacer_cut:" + str(sp1_info["cut"]) + "-" + str(sp1_info["cut"] + 1) + ":0"

        # Beacon amplicon, whole beacon insertion, w/ flank 10bp
        beacon_qw1 = "Beacon:beacon_whole:" + str(sp1_info["cut"] + 1) + "-" + str(
            sp1_info["cut"] + len(beacon)) + ":10"

        subprocess.call(
            "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT,Beacon --guide_seq %s --name %s --output_folder %s "
            "--min_frequency_alleles_around_cut_to_plot 0.05 --write_detailed_allele_table --needleman_wunsch_gap_extend 0 "
            "--trim_sequences --trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:10:0:true "
            "--place_report_in_output_folder --n_processes %s --bam_output --suppress_report %s " % (
                r1, r2, wt_amplicon + "," + beacon_amplicon, sp1_info["seq"], name, output, ncpu, cs2),
            stderr=job_fh, stdout=job_fh, shell=True)

        cs2_stats.update(
            window_quantification(os.path.join(output, "CRISPResso_on_" + name), [wt_qw1, beacon_qw1]))

    pd.Series(cs2_stats).to_json(os.path.join(output, "CRISPResso_on_" + name, "CRISPResso_benchling_stats.json"))
    cs2_stats.update(aaanid=aaan_id, ppid=pp_id)
    pd.Series(cs2_stats).to_json(os.path.join(output, "CRISPResso_on_" + name, "CRISPResso_quilt_stats.json"))

    # plot
    if sp1_info:
        subprocess.call(
            "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a WT --min_freq 0.01 "
            "--plot_center %s --plot_left %s --plot_right %s --plot_cut_point" % (
                os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), sp1_info["cut"] - 1,
                sp1_info["cut"], len(wt_amplicon) - sp1_info["cut"]), stderr=job_fh, stdout=job_fh, shell=True)

        subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s -n 10000000000" % (
            os.path.join(output, "CRISPResso_on_" + name), "WT", wt_qw1), stderr=job_fh, stdout=job_fh, shell=True)

        subprocess.call(
            "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a Beacon --min_freq 0.01 "
            "--plot_center %s --plot_left %s --plot_right %s --plot_cut_point" % (
                os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), sp1_info["cut"] - 1,
                sp1_info["cut"],
                len(beacon_amplicon) - sp1_info["cut"]), stderr=job_fh, stdout=job_fh, shell=True)

        subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s -n 10000000000" % (
            os.path.join(output, "CRISPResso_on_" + name), "Beacon", beacon_qw1), stderr=job_fh, stdout=job_fh, shell=True)

    else:
        subprocess.call(
            "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a WT --min_freq 0.01 "
            "--plot_center %s --plot_left %s --plot_right %s --plot_cut_point" % (
                os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), 0, 1, len(wt_amplicon) - 1),
            stderr=job_fh, stdout=job_fh, shell=True)

        subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -n 10000000000" % (
            os.path.join(output, "CRISPResso_on_" + name), "WT"), stderr=job_fh, stdout=job_fh, shell=True)

    job_fh.close()
    amplicon_fh.close()


if __name__ == "__main__":
    main()
