#!/usr/bin/env python
"""
On-Target analysis for Amp-seq
Information from Benchling
Input AA/SG ID
Input WT corridnates
"""

import argparse
import glob
from utils.base import *
from utils.common_functions import *


def main():
    # Parse command line options
    parser = argparse.ArgumentParser(description='Process meta data', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", help='Path to fastq')
    parser.add_argument("-g1", help='spacer 1')
    parser.add_argument("-g2", help='spacer 2')
    parser.add_argument("-a", help='WT corridnates, eg. hg38:chr1:1-200:+')
    parser.add_argument("-s", help='specifiy a sample in the fastq path', default=None)
    parser.add_argument("-p", help='Number of CPUs to use', default=4)
    parser.add_argument("-o", help='Output folder', default="./")
    parser.add_argument("-cs2", help='CRISPRESSO2 parameters', default="")
    parser.add_argument("-beacon", help='overwrite the beacon sequences', default="")

    args = parser.parse_args()
    fastq = args.i
    name = args.s
    sp1_seq = args.g1
    sp2_seq = args.g2
    ncpu = int(args.p)
    cs2 = args.cs2
    output = args.o
    beacon = args.beacon

    # create output folder
    os.makedirs(os.path.join(output, "cs2_alignment_html"), exist_ok=True)

    # amplicon information file
    amplicon_fh = open(os.path.join(output, name + ".amplicon.txt"), 'w')

    genome_build, chr, region, target_strand = args.a.split(":")
    wt_start, wt_end = region.split("-")

    cs2_stats = {"aaanid": args.g1, "ppid": args.a, "samplename": name}

    # skip if no sample name or no fastq
    if not name or len(glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R1_*")) == 0:
        print("Sample is not in the fastq path")

    # reference genome
    genome_build = re.sub(".*/", "", genome_build)
    genome_fa = "/home/ubuntu/annotation/2bit/" + genome_build + ".2bit"

    # get r1 and r2 fastq
    if target_strand == "antisense" or target_strand == "-":
        r1 = glob.glob(os.path.abspath(fastq) + "/" + name + "*_R2_*")[0]
        r2 = glob.glob(os.path.abspath(fastq) + "/" + name + "*_R1_*")[0]
    else:
        r1 = glob.glob(os.path.abspath(fastq) + "/" + name + "*_R1_*")[0]
        r2 = glob.glob(os.path.abspath(fastq) + "/" + name + "*_R2_*")[0]

    # sample job log
    job_fh = open(os.path.join(output, name + ".job.log"), 'wb')

    # WT amplicon
    wt_amplicon = get_seq(genome_fa, chr, int(wt_start), int(wt_end), target_strand)
    amplicon_fh.write(name + "\tWT\t" + wt_amplicon + "\n")

    sp1_info = {}

    # Beacon amplicon
    sp1_info = get_cut_site(wt_amplicon, sp1_seq)
    sp2_info = get_cut_site(wt_amplicon, sp2_seq)
    if sp1_info["cut"] > sp2_info["cut"]:
        sp1_info = get_cut_site(wt_amplicon, sp2_seq)
        sp2_info = get_cut_site(wt_amplicon, sp1_seq)

    # beacon seq
    beacon_amplicon = wt_amplicon[0:sp1_info["cut"]] + beacon + wt_amplicon[sp2_info["cut"]:]
    amplicon_fh.write(name + "\tBeacon\t" + beacon_amplicon + "\n")

    # define quantification window
    # WT amplicon, spacer cutting 2bp
    wt_qw1 = "WT:spacer1_cut:" + str(sp1_info["cut"]) + "-" + str(sp1_info["cut"] + 1) + ":0"
    # WT amplicon, spacer2 cutting 2bp
    wt_qw2 = "WT:spacer2_cut:" + str(sp2_info["cut"]) + "-" + str(sp2_info["cut"] + 1) + ":0"
    # Beacon amplicon, whole beacon insertion, w/ flank 10bp
    beacon_qw1 = "Beacon:beacon_whole:" + str(sp1_info["cut"] + 1) + "-" + str(sp1_info["cut"] + len(beacon)) + ":10"

    subprocess.call(
        "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT,Beacon --guide_seq %s --name %s --output_folder %s "
        "--min_frequency_alleles_around_cut_to_plot 0.05 --write_detailed_allele_table --needleman_wunsch_gap_extend 0 "
        "--trim_sequences  --trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:10:0:true "
        "--place_report_in_output_folder --n_processes %s --bam_output --suppress_report %s " % (
            r1, r2, wt_amplicon + "," + beacon_amplicon, sp1_info["seq"] + "," + sp2_info["seq"], name, output, ncpu, cs2),
        stderr=job_fh, stdout=job_fh, shell=True)

    cs2_stats.update(window_quantification(os.path.join(output, "CRISPResso_on_" + name), [wt_qw1, wt_qw2, beacon_qw1]))

    pd.Series(cs2_stats).to_json(os.path.join(output, "CRISPResso_on_" + name, "CRISPResso_benchling_stats.json"))
    cs2_stats.update(aaanid=args.g1, ppid=args.a)
    pd.Series(cs2_stats).to_json(os.path.join(output, "CRISPResso_on_" + name, "CRISPResso_quilt_stats.json"))

    # plot
    subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s -b %s -n 10000000000" % (
        os.path.join(output, "CRISPResso_on_" + name), "WT", wt_qw1, wt_qw2), stderr=job_fh, stdout=job_fh, shell=True)
    subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s -n 10000000000" % (
            os.path.join(output, "CRISPResso_on_" + name), "Beacon", beacon_qw1), stderr=job_fh, stdout=job_fh, shell=True)

    job_fh.close()
    amplicon_fh.close()


if __name__ == "__main__":
    main()
