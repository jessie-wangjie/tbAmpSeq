#!/usr/bin/env python
"""
Amp-Seq for UMASS settings, single atgRNA with PBS2
Information from Excel sheet
"""

import argparse
import glob
import datetime
import os.path
from logging import warning
from utils.base import *
from utils.common_functions import *


def main():
    # Parse command line options
    parser = argparse.ArgumentParser(description='Process meta data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", help='metasheet')
    parser.add_argument("-i", help='Path to fastq')
    parser.add_argument("-s", help='specific sample to analysis', default=None)
    parser.add_argument("-p", help='Number of CPUs to use', default=4)
    parser.add_argument("-o", help='Output folder', default="./")
    parser.add_argument("-cs2", help='CRISPRESSO2 parameters', default="")

    args = parser.parse_args()
    fastq = args.i
    sample = args.s
    ncpu = int(args.p)
    cs2 = args.cs2
    output = args.o

    tbid = os.path.basename(os.path.abspath(fastq))
    print(tbid)
    # read in metasheet
    info = pd.read_excel(args.m)

    try:
        os.makedirs(output)
    except:
        warning('Folder %s already exists.' % output)

    amplicon_fh = open(os.path.join(output, tbid + ".amplicon.txt"), 'w')
    project_fh = open(os.path.join(output, tbid + ".status.txt"), 'w')

    run_start = str(datetime.datetime.now())
    if os.path.exists(os.path.join(output, tbid + ".run.json")):
        ngs_stats = pd.read_json(os.path.join(output, tbid + ".run.json"), typ="series")
    else:
        ngs_stats = pd.Series(dtype="object")

    for i in info.index:
        s = info.loc[i]

        name = s["miseq_sample_name"]
        aaan_id = s["aaanid"]
        pp_id = s["ppid"]

        cs2_stats = {"aaanid": s["aaanid"], "ppid": s["ppid"], "samplename": s["samplename"],
                     "mrna_batch_id": s["mrna_batch_id"], "modatg_batch_id": s["modatg_batch_id"],
                     "primary_cell_lot_id": s["primary_cell_lot_id"], "lnp_batch_id": s["lnp_batch_id"],
                     "plate": s["plate"], "well": s["well"], "miseq_sample_name": name,
                     "genomics_ampseq_project_queue": tbid}
        print([name, aaan_id, pp_id])

        if pd.isna(name):
            continue

        if not name or len(glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R1_*")) == 0:
            continue

        if sample and name != sample:
            continue

        # Get primer information
        # special case PRIP874
        cur.execute(
            "select p1.chromosome, p2.start, p1.end, p1.genome_build, target_gene.direction_of_transcription from primer_pair "
            "join primer as p1 on p1.id = primer_pair.forward_primer "
            "join primer as p2 on p2.id = primer_pair.reverse_primer "
            "join target_gene on target_gene.id = p1.gene_or_target_name "
            "where primer_pair.file_registry_id$ = %s", [pp_id])
        target_chr, wt_start, wt_end, genome_build, target_strand = cur.fetchone()

        # reference genome
        genome_build = re.sub(".*/", "", genome_build)
        genome_fa = "/home/ubuntu/annotation/2bit/" + genome_build + ".2bit"

        # get r1 and r2 fastq
        # if target_strand == "antisense" or target_strand == "-":
        #    r1 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R2_*")[0]
        #     r2 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R1_*")[0]
        # else:
        r1 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R1_*")[0]
        r2 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R2_*")[0]

        job_fh = open(os.path.join(output, name + ".job.log"), 'wb')

        # WT amplicon
        wt_amplicon = get_seq(genome_fa, target_chr, wt_start, wt_end, target_strand)
        amplicon_fh.write(name + "\tWT\t" + wt_amplicon + "\n")

        sp1_info = {}

        # Beacon amplicon
        # atgRNA-ngRNA
        if not aaan_id:
            subprocess.call(
                "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s "
                "--needleman_wunsch_gap_extend 0 %s --bam_output --suppress_report --trim_sequences "
                "--trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:7:0:true" % (
                    r1, r2, wt_amplicon, name, output, ncpu, cs2), stderr=job_fh, stdout=job_fh, shell=True)

        elif aaan_id.startswith("AN"):
            # get spacer sequences, beacon sequences, ngRNA sequences, pbs2 seqeunces
            sp1_info = get_cut_site(wt_amplicon, s["sp_seq"])
            ng_info = get_cut_site(wt_amplicon, s["ng_seq"])

            pbs2_info = get_cut_site(wt_amplicon, s["pbs2_seq"])
            beacon = get_beacon_seq(s["beacon_seq"], sp1_info["strand"])

            # beacon seq
            beacon_amplicon = wt_amplicon[0:sp1_info["cut"]] + beacon + wt_amplicon[pbs2_info["3P"] - 1:]
            amplicon_fh.write(name + "\tBeacon\t" + beacon_amplicon + "\n")

            # define quantification window
            # WT amplicon, spacer cutting 2bp
            wt_qw1 = "WT:spacer_cut:" + str(sp1_info["cut"]) + "-" + str(sp1_info["cut"] + 1) + ":0"

            # WT amplicon, ngRNA cutting 2bp
            wt_qw2 = "WT:ng_cut:" + str(ng_info["cut"]) + "-" + str(ng_info["cut"] + 1) + ":0"

            # Beacon amplicon, whole beacon insertion, w/ flank 10bp
            beacon_qw1 = "Beacon:beacon_whole:" + str(sp1_info["cut"] + 1) + "-" + str(
                sp1_info["cut"] + len(beacon)) + ":10"

            # Beacon amplicon, RT 5'
            if pbs2_info["strand"] == "+":
                beacon_qw2 = "Beacon:RT_5P:" + str(pbs2_info["5P"] - 1) + "-" + str(pbs2_info["5P"]) + ":0"
            else:
                beacon_qw2 = "Beacon:RT_5P:" + str(pbs2_info["5P"]) + "-" + str(pbs2_info["5P"] + 1) + ":0"

            # Beacon amplicon, ngRNA cutting 2bp
            ng_info = get_cut_site(beacon_amplicon, s["ng_seq"])
            beacon_qw3 = "Beacon:ng_cut:" + str(ng_info["cut"]) + "-" + str(ng_info["cut"] + 1) + ":0"

            subprocess.call(
                "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT,Beacon --guide_seq %s "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s "
                "--needleman_wunsch_gap_extend 0 %s --bam_output --suppress_report --trim_sequences "
                "--trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:7:0:true" % (
                    r1, r2, wt_amplicon + "," + beacon_amplicon, sp1_info["seq"] + "," + ng_info["seq"], name, output,
                    ncpu, cs2), stderr=job_fh, stdout=job_fh, shell=True)

            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s -qw %s -qw %s" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    wt_qw1, wt_qw2, beacon_qw1, beacon_qw2, beacon_qw3), stderr=job_fh, stdout=job_fh, shell=True)

            cs2_stats.update(window_quantification(os.path.join(output, "CRISPResso_on_" + name),
                                                   [wt_qw1, wt_qw2, beacon_qw1, beacon_qw2, beacon_qw3]))

        pd.concat([pd.Series(cs2_stats), ngs_stats]).to_json(
            os.path.join(output, "CRISPResso_on_" + name, "CRISPResso_stats.json"))

        # write sample status
        with open(os.path.join(output, "CRISPResso_on_" + name, "CRISPResso_status.txt"), 'r') as f:
            status = f.readline()
            for line in f:
                if "ERROR" in line:
                    status = line
                    p = subprocess.Popen(
                        "samtools view %s | cut -f10 | sort | uniq -c | sort -n | tail -1" % (
                            os.path.join(output, "CRISPResso_on_" + name, "CRISPResso_output.bam")), stdout=subprocess.PIPE, shell=True)
                    status = status.rstrip() + "\t" + p.communicate()[0].decode('utf-8').lstrip()
            project_fh.write(name + "\t" + status)

        # plot
        if sp1_info:
            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a WT --plot_center %s "
                "--plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    sp1_info["cut"] - 1, sp1_info["cut"], len(wt_amplicon) - sp1_info["cut"]), stderr=job_fh,
                stdout=job_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "WT", wt_qw1), stderr=job_fh, stdout=job_fh,
                            shell=True)
        else:
            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a WT --plot_center %s "
                "--plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), 0, 1,
                    len(wt_amplicon) - 1), stderr=job_fh, stdout=job_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "WT"), stderr=job_fh, stdout=job_fh, shell=True)

        if aaan_id and aaan_id.startswith("PN"):
            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a Prime-edited --plot_center %s "
                "--plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    sp1_info["cut"] - 1, sp1_info["cut"], len(beacon_amplicon) - sp1_info["cut"]), stderr=job_fh,
                stdout=job_fh, shell=True)
            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a Scaffold-incorporated --plot_center %s "
                "--plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    sp1_info["cut"] - 1, sp1_info["cut"], len(beacon_amplicon) - sp1_info["cut"]), stderr=job_fh,
                stdout=job_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "Prime-edited", beacon_qw1), stderr=job_fh,
                            stdout=job_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "Scaffold-incorporated", beacon_qw1), stderr=job_fh,
                            stdout=job_fh, shell=True)

        elif aaan_id and (not aaan_id.startswith("SG")):
            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a Beacon --plot_center %s "
                "--plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    sp1_info["cut"] - 1, sp1_info["cut"], len(beacon_amplicon) - sp1_info["cut"]), stderr=job_fh,
                stdout=job_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "Beacon", beacon_qw1), stderr=job_fh, stdout=job_fh,
                            shell=True)

        job_fh.close()
    amplicon_fh.close()
    project_fh.close()
    ngs_stats["run start"] = run_start
    ngs_stats["run end"] = str(datetime.datetime.now())
    pd.Series(ngs_stats).to_json(os.path.join(output, tbid + ".run.json"))


if __name__ == "__main__":
    main()
