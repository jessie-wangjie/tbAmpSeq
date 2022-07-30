#!/usr/bin/env python
"""
On-Target analysis for Amp-seq and 3Primer-Seq
Information from Benchling
"""

import pandas as pd
import argparse
import subprocess
from utils.common_functions import *
from utils.base import *
import os
import glob
import re


def main():
    # Parse command line options
    parser = argparse.ArgumentParser(description='Process meta data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", help='Path to the meta table')
    parser.add_argument("-i", help='Path to fastq')
    parser.add_argument("-p", help='Number of CPUs to use', default=4)
    parser.add_argument("-o", help='Output folder')

    args = parser.parse_args()
    fastq = args.i
    ncpu = int(args.p)
    info = pd.read_excel(args.m, skiprows=21)
    output = args.o

    amplicon_fh = open(os.path.join(output, os.path.basename(fastq) + ".amplicon.txt"), 'w')

    # Read from metasheet
    for i in info.index:
        sample = info.loc[i]

        # Finish processing all the samples
        if sample["Sample name"] == "PROTOCOLS":
            break

        # Skip empty lines
        if pd.isna(sample["Sample name"]):
            continue

        # Get primer information
        if pd.isna(sample["PP ID"]):
            cur.execute(
                "select primer.chromosome, primer.genome_build, target_gene.direction_of_transcription from dna_oligo "
                "join primer on primer.id=dna_oligo.id "
                "join target_gene on target_gene.id = primer.gene_or_target_name "
                "where dna_oligo.bases = %s", [sample["Forward Primer Sequence"]])
            target_chr, genome_build, target_strand = cur.fetchone()
            print(target_chr)
            print(target_strand)
            reference_index = "/home/ubuntu/annotation/bwa_index/" + genome_build
            fp_info = align_primer(sample["Forward Primer Sequence"], reference_index, target_chr, "CACTCTTTCCCTACACGACGCTCTTCCGATCT")
            rp_info = align_primer(sample["Reverse Primer Sequence"], reference_index, target_chr, "GGAGTTCAGACGTGTGCTCTTCCGATCT")

            print(fp_info)
            wt_start = fp_info["start"]
            wt_end = rp_info["end"]
        else:
            cur.execute(
                "select p1.chromosome, p1.start, p2.end, p1.genome_build, target_gene.direction_of_transcription from primer_pair "
                "join primer as p1 on p1.id = primer_pair.forward_primer "
                "join primer as p2 on p2.id = primer_pair.reverse_primer "
                "join target_gene on target_gene.id = p1.gene_or_target_name "
                "where primer_pair.file_registry_id$ = %s", [sample["PP ID"]])
            target_chr, wt_start, wt_end, genome_build, target_strand = cur.fetchone()

        beacon_amplicon = ""
        cargo_amplicon = ""
        assay = "AmpSeq"
        name = sample["MiSeq samplesheet name (no space, no underscore)"]

        # reference genome
        genome_fa = "/home/ubuntu/annotation/2bit/" + genome_build + ".2bit"
        cargo_fa = "/home/ubuntu/annotation/2bit/PL224.nanoluc.2bit"

        # get r1 and r2 fastq
        print(name)
        r1 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R1_*")[0]
        r2 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R2_*")[0]

        # WT amplicon
        wt_amplicon = get_seq(genome_fa, target_chr, wt_start, wt_end, target_strand)
        amplicon_fh.write(name + "\tWT\t" + wt_amplicon + "\n")

        sp1_info = {}
        rt_info = {}
        wt_qw1 = ""

        # Beacon amplicon
        # locate atgRNA, ngRNA in the amplicon
        if sample["atgRNA pairing type"] == "single":
            # get spacer sequences, beacon sequences, ngRNA sequences
            cur.execute("select sp.bases, beacon.bases, ng.bases, atgrna.rt_coordinate, atg.bases from atg_ng "
                        "join modified_rna as m1 on m1.id=atg_ng.atgrna "
                        "join modified_rna as m2 on m2.id=atg_ng.ngrna "
                        "join atgrna on atgrna.id=m1.rna "
                        "join ngrna on ngrna.id=m2.rna "
                        "join dna_oligo as sp on sp.id=atgrna.spacer "
                        "join dna_sequence as beacon on beacon.id=atgrna.beacon "
                        "join dna_oligo as ng on ng.id=ngrna.spacer "
                        "join dna_sequence as atg on atg.id=atgrna.id "
                        "where atg_ng.file_registry_id$ = %s", [sample["AA/AN ID"]])
            sp_seq, beacon_seq, ng_seq, rt_coord, atg_seq = cur.fetchone()
            sp1_info = get_cut_site(wt_amplicon, sp_seq)
            ng_info = get_cut_site(wt_amplicon, ng_seq)
            coord = re.match(r"\[(.*), (.*)\]", rt_coord)
            print(coord)
            print(coord.group(1))
            print(coord.group(2))
            rt_info = get_cut_site(wt_amplicon, atg_seq[int(coord.group(1))-1:int(coord.group(2))])
            print(rt_info["seq"])
            beacon = get_beacon_seq(beacon_seq, sp1_info["strand"])

            # beacon seq
            beacon_amplicon = wt_amplicon[0:sp1_info["cut"]] + beacon + wt_amplicon[sp1_info["cut"]:]
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
            if rt_info["strand"] == "+":
                beacon_qw2 = "Beacon:RT_5P:" + str(rt_info["5P"] - 1) + "-" + str(rt_info["5P"]) + ":0"
            else:
                beacon_qw2 = "Beacon:RT_5P:" + str(rt_info["5P"]) + "-" + str(rt_info["5P"] + 1) + ":0"

            # Beacon amplicon, ngRNA cutting 2bp
            ng_info = get_cut_site(beacon_amplicon, ng_seq)
            beacon_qw3 = "Beacon:ng_cut:" + str(ng_info["cut"]) + "-" + str(ng_info["cut"] + 1) + ":0"

        # locate atgRNA, atgRNA in the amplicon
        elif sample["atgRNA pairing type"] == "dual":
            # Get spacers information
            cur.execute( "select sp1.bases, sp2.bases, beacon1.bases, beacon2.bases from atg_atg "
                         "join modified_rna as m1 on m1.id = atg_atg.atg1 "
                         "join modified_rna as m2 on m2.id = atg_atg.atg2 "
                         "join atgrna as a1 on a1.id = m1.rna "
                         "join atgrna as a2 on a2.id = m2.rna "
                         "join dna_oligo as sp1 on sp1.id = a1.spacer "
                         "join dna_oligo as sp2 on sp2.id = a2.spacer "
                         "join dna_sequence as beacon1 on beacon1.id = a1.beacon "
                         "join dna_sequence as beacon2 on beacon2.id = a2.beacon "
                         "where atg_atg.file_registry_id$ = %s", [sample["AA/AN ID"]])
            sp1_seq, sp2_seq, beacon1_seq, beacon2_seq = cur.fetchone()

            sp1_info = get_cut_site(wt_amplicon, sp1_seq)
            sp2_info = get_cut_site(wt_amplicon, sp2_seq)

            beacon = get_beacon_seq(beacon1_seq, sp1_info["strand"], beacon2_seq, sp2_info["strand"])

            # beacon seq
            beacon_amplicon = wt_amplicon[0:sp1_info["cut"]] + beacon + wt_amplicon[sp2_info["cut"]:]
            amplicon_fh.write(name + "\tBeacon\t" + beacon_amplicon + "\n")

            # define quantification window
            # WT amplicon, spacer cutting 2bp
            wt_qw1 = "WT:spacer1_cut:" + str(sp1_info["cut"]) + "-" + str(sp1_info["cut"] + 1) + ":0"
            # WT amplicon, spacer2 cutting 2bp
            wt_qw2 = "WT:spacer2_cut:" + str(sp2_info["cut"]) + "-" + str(sp2_info["cut"] + 1) + ":0"
            # Beacon amplicon, whole beacon insertion, w/ flank 10bp
            beacon_qw1 = "Beacon:beacon_whole:" + str(sp1_info["cut"] + 1) + "-" + str(
                sp1_info["cut"] + len(beacon)) + ":10"

        else:
            assay_type = "Control"

        # cargo seq
        if assay == "3P":
            beacon_info = align_primer(beacon[:18], "/home/ubuntu/annotation/bwa_index/PL224.nanoluc")
            cargo_amplicon = wt_amplicon[0:sp1_info["cut"]] + get_seq(cargo_fa, rp2_info["chr"],
                                                                         beacon_info["start"], rp2_info["end"], "+")
            amplicon_fh.write(name + "\tCargo\t" + cargo_amplicon + "\n")

            # define quantification window
            cargo_qw1 = "Cargo:spacer_cut:" + str(sp1_info["cut"]) + "-" + str(sp1_info["cut"] + 1) + ":0"
            cargo_qw2 = "Cargo:AttL1:" + str(sp1_info["cut"] + 1) + "-" + str(sp1_info["cut"] + 30) + ":0"

        # Run Crispresso2
        error_fh = open(os.path.join(output, name + ".job.log"), 'wb')

        ## mkdir folder
        preprocess_output = os.path.join(output, name + "_preprocess_output")
        if not os.path.exists(preprocess_output):
            os.mkdir(preprocess_output)

        ## merge r1 and r2
        if target_strand == "antisense":
            subprocess.call(
                "flash %s %s --min-overlap 10 --max-overlap 100 --allow-outies -z -d %s" % (r2, r1, preprocess_output),
                shell=True, stderr=error_fh, stdout=error_fh)
        else:
            subprocess.call(
                "flash %s %s --min-overlap 10 --max-overlap 100 --allow-outies -z -d %s" % (r1, r2, preprocess_output),
                shell=True, stderr=error_fh, stdout=error_fh)

        # align to plasmid
        if "Payload ID" in sample and pd.notna(sample["Payload ID"]):
            plasmid_id = {"Pdy0186": 140, "PL249": 140}
            input_file = os.path.join(preprocess_output, "out.extendedFrags.fastq.gz")
            plasmid_sam = os.path.join(preprocess_output, "out.extendedFrags.plasmid.sam")
            unmapped_id = os.path.join(preprocess_output, "out.extendedFrags.unmapped_plasmid.id")
            unmapped_fastq = os.path.join(preprocess_output, "out.extendedFrags.unmapped_plasmid.fastq")
            subprocess.call("bowtie2 --local -p 8 -x %s -U %s -S %s" % (
                "/home/ubuntu/annotation/bowtie2_index/" + sample["Payload ID"], input_file, plasmid_sam), shell=True,
                            stderr=error_fh, stdout=error_fh)
            subprocess.call("samtools view -F4 %s | awk '{ if($4>=%s) print $1}' > %s" % (
                plasmid_sam, plasmid_id[sample["Payload ID"]], unmapped_id), shell=True, stderr=error_fh,
                            stdout=error_fh)
            subprocess.call("samtools view -f4 %s | cut -f1 >> %s" % (plasmid_sam, unmapped_id), shell=True,
                            stderr=error_fh, stdout=error_fh)
            subprocess.call(
                "samtools view -h -N %s %s | samtools fastq -0 %s -" % (unmapped_id, plasmid_sam, unmapped_fastq),
                shell=True, stderr=error_fh, stdout=error_fh)
        else:
            unmapped_fastq = os.path.join(preprocess_output, "out.extendedFrags.fastq.gz")

        # cs2, quantification, plot
        if assay == "Control":
            if wt_qw1:
                subprocess.call(
                    "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name WT --guide_seq %s "
                    "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                    "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s" % (
                        unmapped_fastq, wt_amplicon, ng_info["seq"], name, output, ncpu), stderr=error_fh,
                    stdout=error_fh,
                    shell=True)
                subprocess.call("python utils/parse_quantification_windows.py -f %s -o %s -qw %s" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    wt_qw1),
                                stderr=error_fh, stdout=error_fh, shell=True)
            else:
                subprocess.call(
                    "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name WT "
                    "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                    "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s" % (
                        unmapped_fastq, wt_amplicon, name, output, ncpu), stderr=error_fh, stdout=error_fh, shell=True)

            subprocess.call(
                "python utils/plotCustomAllelePlot.py -f %s -o %s --plot_center %s --plot_left %s "
                "--plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), 0, 1,
                    len(wt_amplicon) - 1), stderr=error_fh, stdout=error_fh, shell=True)
            subprocess.call("allele2html.py -f %s -r %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "WT"), stderr=error_fh, stdout=error_fh, shell=True)

        elif assay == "AmpSeq":
            subprocess.call(
                "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name WT,Beacon --guide_seq %s "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s "
                "--needleman_wunsch_gap_extend 0" % (
                    unmapped_fastq, wt_amplicon + "," + beacon_amplicon, sp1_info["seq"], name, output, ncpu),
                stderr=error_fh, stdout=error_fh, shell=True)
            if sample["atgRNA pairing type"] == "single":
                subprocess.call(
                    "python utils/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s -qw %s" % (
                        os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                        wt_qw1, wt_qw2, beacon_qw1, beacon_qw3), stderr=error_fh, stdout=error_fh, shell=True)
            elif sample["atgRNA pairing type"] == "dual":
                subprocess.call(
                    "python utils/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s" % (
                        os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                        wt_qw1, wt_qw2, beacon_qw1), stderr=error_fh, stdout=error_fh, shell=True)
            subprocess.call(
                "python utils/plotCustomAllelePlot.py -f %s -o %s -a WT --plot_center %s --plot_left %s "
                "--plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    sp1_info["cut"] - 1, sp1_info["cut"], len(wt_amplicon) - sp1_info["cut"]), stderr=error_fh,
                stdout=error_fh, shell=True)
            subprocess.call(
                "python utils/plotCustomAllelePlot.py -f %s -o %s -a Beacon --plot_center %s --plot_left "
                "%s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    sp1_info["cut"] - 1, sp1_info["cut"], len(beacon_amplicon) - sp1_info["cut"]),
                stderr=error_fh,
                stdout=error_fh, shell=True)
            subprocess.call("python utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "WT", wt_qw1), stderr=error_fh, stdout=error_fh,
                            shell=True)
            subprocess.call("python utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "Beacon", beacon_qw1), stderr=error_fh, stdout=error_fh,
                            shell=True)

        elif assay == "3P":
            ## sort reads by RP1 and RP2, with max 2 mismatches
            RP1_fastq = os.path.join(preprocess_output, "RP1.fastq.gz")
            RP2_fastq = os.path.join(preprocess_output, "RP2.fastq.gz")
            noRP_fastq = os.path.join(preprocess_output, "noRP.fastq.gz")
            if gene_strand == "-":
                fq = subprocess.Popen(
                    "/home/ubuntu/software/miniconda3/bin/cutadapt -m 10 -O 10 -e 2.5 -a %s --action=none -o %s "
                    "--untrimmed-output - %s | /home/ubuntu/software/miniconda3/bin/cutadapt -m 10 -O 10 -e 2.5 -a %s "
                    "--action=none --untrimmed-output %s -o %s -" % (
                        reverse_complement(FP_info["seq"]), RP1_fastq, unmapped_fastq,
                        reverse_complement(RP2_info["seq"]),
                        noRP_fastq, RP2_fastq), shell=True, stdout=error_fh, stderr=error_fh)
            else:
                fq = subprocess.Popen(
                    "/home/ubuntu/software/miniconda3/bin/cutadapt -m 10 -O 10 -e 2.5 -a %s --action=none -o %s "
                    "--untrimmed-output - %s | /home/ubuntu/software/miniconda3/bin/cutadapt -m 10 -O 10 -e 2.5 -a %s "
                    "--action=none --untrimmed-output %s -o %s -" % (
                        reverse_complement(RP1_info["seq"]), RP1_fastq, unmapped_fastq,
                        reverse_complement(RP2_info["seq"]),
                        noRP_fastq, RP2_fastq), shell=True, stdout=error_fh, stderr=error_fh)

            subprocess.call(
                "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name WT,Beacon --guide_seq %s "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s" % (
                    RP1_fastq, wt_amplicon + "," + beacon_amplicon, sp1_info["seq"], name + "_WT_Beacon", output,
                    ncpu),
                stderr=error_fh, stdout=error_fh, shell=True)
            subprocess.call(
                "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name Cargo "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s" % (
                    RP2_fastq, cargo_amplicon, name + "_Cargo", output, ncpu), stderr=error_fh, stdout=error_fh,
                shell=True)

            if assaytype == "single_atg":
                subprocess.call(
                    "python utils/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s -qw %s "
                    "-qw %s" % (
                        os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"),
                        os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"), wt_qw1, wt_qw2, beacon_qw1,
                        beacon_qw2, beacon_qw3), stderr=error_fh, stdout=error_fh, shell=True)
            elif assaytype == "dual_atg":
                subprocess.call(
                    "python utils/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s" % (
                        os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"),
                        os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"), wt_qw1, wt_qw2, beacon_qw1),
                    stderr=error_fh, stdout=error_fh, shell=True)

            subprocess.call("python utils/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s" % (
                os.path.join(output, "CRISPResso_on_" + name + "_Cargo"),
                os.path.join(output, "CRISPResso_on_" + name + "_Cargo"), cargo_qw1, cargo_qw2), stderr=error_fh,
                            stdout=error_fh, shell=True)

            subprocess.call(
                "python utils/plotCustomAllelePlot.py -f %s -o %s -a WT --plot_center %s --plot_left %s "
                "--plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"),
                    os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"), sp1_info["cut"] - 1,
                    sp1_info["cut"], len(wt_amplicon) - sp1_info["cut"]), stderr=error_fh, stdout=error_fh,
                shell=True)
            subprocess.call(
                "python utils/plotCustomAllelePlot.py -f %s -o %s -a Beacon --plot_center %s --plot_left "
                "%s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"),
                    os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"), sp1_info["cut"] - 1,
                    sp1_info["cut"], len(beacon_amplicon) - sp1_info["cut"]), stderr=error_fh, stdout=error_fh,
                shell=True)
            subprocess.call(
                "python utils/plotCustomAllelePlot.py -f %s -o %s -a Cargo --plot_center %s --plot_left %s "
                "--plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name + "_Cargo"),
                    os.path.join(output, "CRISPResso_on_" + name + "_Cargo"), sp1_info["cut"] - 1,
                    sp1_info["cut"],
                    len(cargo_amplicon) - sp1_info["cut"]), stderr=error_fh, stdout=error_fh, shell=True)
            subprocess.call("python utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"), "WT", wt_qw1), stderr=error_fh,
                            stdout=error_fh, shell=True)
            subprocess.call("python utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"), "Beacon", beacon_qw1), stderr=error_fh,
                            stdout=error_fh, shell=True)
            subprocess.call("python utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name + "_Cargo"), "Cargo", cargo_qw2), stderr=error_fh,
                            stdout=error_fh, shell=True)
        error_fh.close()
    amplicon_fh.close()


if __name__ == "__main__":
    main()
