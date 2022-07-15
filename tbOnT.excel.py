#!/usr/bin/env python
"""
Tome 3P/AmpSeq meta data
"""

import pandas as pd
import argparse
import subprocess
import re
import os
import glob

def main():
    # Parse command line options
    parser = argparse.ArgumentParser(description='Process meta data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", help='Path to the meta table')
    parser.add_argument("-i", help='Path to fastq')
    parser.add_argument("-p", help='Number of CPUs to use', default=4)
    parser.add_argument("-o", help='Output folder', default=4)

    args = parser.parse_args()
    fastq = args.i
    ncpu = int(args.p)
    info = pd.read_excel(args.m, skiprows=21)
    output = args.o

    amplicon_fh = open(os.path.join(output, os.path.basename(fastq) + ".amplicon.txt"), 'w')

    illumina = {"P5": "CACTCTTTCCCTACACGACGCTCTTCCGATCT", "P7": "GGAGTTCAGACGTGTGCTCTTCCGATCT"}

    for i in info.index:
        sample = info.loc[i]

        if sample["Sample name"] == "PROTOCOLS":
            break

        if pd.isna(sample["Sample name"]):
            continue

        wt_amplicon = ""
        beacon_amplicon = ""
        cargo_amplicon = ""
        assay = "AmpSeq"
        gene_strand = sample["direction of transcription of target gene"]
        name = sample["MiSeq samplesheet name (no space, no underscore)"]

        # reference genome
        reference_index = "/home/ubuntu/annotation/bwa_index/" + re.sub(".* ", "", sample["Genome Build"])
        genome_fa = "/home/ubuntu/annotation/2bit/" + re.sub(".* ", "", sample["Genome Build"]) + ".2bit"
        cargo_fa = "/home/ubuntu/annotation/2bit/PL224.nanoluc.2bit"

        # get R1 and R2
        print(name)
        R1 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R1_*")[0]
        R2 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R2_*")[0]

        # locate the FP and RP in the genome
        if "Forward Primer link to Illumina adapter" not in sample:
            FP_adapter = "P5"
        else:
            FP_adapter = sample["Forward Primer link to Illumina adapter"]
        FP_info = align_primer(sample["Forward Primer Sequence"], reference_index, illumina[FP_adapter])
        if "Reverse Primer Sequence" in sample:
            if "Reverse Primer link to Illumina adapter" not in sample:
                RP1_adapter = "P7"
            else:
                RP1_adapter = sample["Reverse Primer link to Illumina adapter"]
            RP1_info = align_primer(sample["Reverse Primer Sequence"], reference_index, illumina[RP1_adapter])
        else:
            if "Reverse Primer 1 link to Illumina adapter" not in sample:
                RP1_adapter = "P7"
            else:
                RP1_adapter = sample["Reverse Primer 1 link to Illumina adapter"]
            RP1_info = align_primer(sample["Reverse Primer 1 Sequence"], reference_index, illumina[RP1_adapter])

        if "Reverse Primer 2 Sequence" in sample:
            if sample["characteristics: treatment"] == "Treated":
                assay = "3P"
            if "Reverse Primer 2 link to Illumina Adapater" not in sample:
                RP2_adapter = "P7"
            else:
                RP2_adapter = sample["Reverse Primer 2 link to Illumina Adapater"]
            RP2_info = align_primer(sample["Reverse Primer 2 Sequence"],
                                    "/home/ubuntu/annotation/bwa_index/PL224.nanoluc", illumina[RP2_adapter])

        wt_amplicon = get_seq(genome_fa, FP_info["chr"], FP_info["start"], RP1_info["end"], gene_strand)
        amplicon_fh.write(name + "\tWT\t" + wt_amplicon + "\n")

        spacer_info = {}
        rt_info = {}
        wt_qw1 = ""
        # Beacon amplicon
        # locate atgRNA, ngRNA in the amplicon
        if "atg-ng pair ID" in sample and pd.notna(sample["atg-ng pair ID"]):
            type = "single_atg"
            spacer_info = get_cut_site(wt_amplicon, sample["atgRNA spacer sequence"])
            ng_info = get_cut_site(wt_amplicon, sample["ngRNA spacer sequence"])
            if pd.notna(sample["atgRNA RT seq"]):
                rt_info = get_cut_site(wt_amplicon, sample["atgRNA RT seq"])
            beacon = get_beacon_seq(sample["atgRNA beacon seq"], spacer_info["strand"])
            # beacon seq
            beacon_amplicon = wt_amplicon[0:spacer_info["cut"]] + beacon + wt_amplicon[spacer_info["cut"]:]
            amplicon_fh.write(name + "\tBeacon\t" + beacon_amplicon + "\n")

            # define quantification window
            # WT amplicon, spacer cutting 2bp
            wt_qw1 = "WT:spacer_cut:" + str(spacer_info["cut"]) + "-" + str(spacer_info["cut"] + 1) + ":0"
            # WT amplicon, ngRNA cutting 2bp
            wt_qw2 = "WT:ng_cut:" + str(ng_info["cut"]) + "-" + str(ng_info["cut"] + 1) + ":0"
            # Beacon amplicon, whole beacon insertion, w/ flank 10bp
            beacon_qw1 = "Beacon:beacon_whole:" + str(spacer_info["cut"] + 1) + "-" + str(
                spacer_info["cut"] + len(beacon)) + ":10"
            # Beacon amplicon, RT 5'
            if not rt_info:
                beacon_qw2 = ""
            elif rt_info["strand"] == "+":
                beacon_qw2 = "Beacon:RT_5P:" + str(rt_info["5P"] - 1) + "-" + str(rt_info["5P"]) + ":0"
            else:
                beacon_qw2 = "Beacon:RT_5P:" + str(rt_info["5P"]) + "-" + str(rt_info["5P"] + 1) + ":0"
            # Beacon amplicon, ngRNA cutting 2bp
            ng_info = get_cut_site(beacon_amplicon, sample["ngRNA spacer sequence"])
            beacon_qw3 = "Beacon:ng_cut:" + str(ng_info["cut"]) + "-" + str(ng_info["cut"] + 1) + ":0"

        # locate atgRNA, atgRNA in the amplicon
        elif "atg-atg pair ID" in sample and pd.notna(sample["atg-atg pair ID"]):
            type = "dual_atg"
            spacer_info = get_cut_site(wt_amplicon, sample["atgRNA1 spacer sequence"])
            spacer2_info = get_cut_site(wt_amplicon, sample["atgRNA2 spacer sequence"])
            beacon = get_beacon_seq(sample["atgRNA1 beacon seq"], spacer_info["strand"], sample["atgRNA2 beacon seq"],
                                    spacer2_info["strand"])
            # beacon seq
            beacon_amplicon = wt_amplicon[0:spacer_info["cut"]] + beacon + wt_amplicon[spacer2_info["cut"]:]
            amplicon_fh.write(name + "\tBeacon\t" + beacon_amplicon + "\n")

            # define quantification window
            # WT amplicon, spacer cutting 2bp
            wt_qw1 = "WT:spacer1_cut:" + str(spacer_info["cut"]) + "-" + str(spacer_info["cut"] + 1) + ":0"
            # WT amplicon, spacer2 cutting 2bp
            wt_qw2 = "WT:spacer2_cut:" + str(spacer2_info["cut"]) + "-" + str(spacer2_info["cut"] + 1) + ":0"
            # Beacon amplicon, whole beacon insertion, w/ flank 10bp
            beacon_qw1 = "Beacon:beacon_whole:" + str(spacer_info["cut"] + 1) + "-" + str(
                spacer_info["cut"] + len(beacon)) + ":10"

        else:
            assay = "Control"
            if "ngRNA ID" in sample and pd.notna(sample["ngRNA ID"]):
                ng_info = get_cut_site(wt_amplicon, sample["ngRNA spacer sequence"])
                wt_qw1 = "WT:ng_cut:" + str(ng_info["cut"]) + "-" + str(ng_info["cut"] + 1) + ":0"

        # cargo seq
        if assay == "3P":
            beacon_info = align_primer(beacon[:18], "/home/ubuntu/annotation/bwa_index/PL224.nanoluc")
            cargo_amplicon = wt_amplicon[0:spacer_info["cut"]] + get_seq(cargo_fa, RP2_info["chr"],
                                                                         beacon_info["start"], RP2_info["end"], "+")
            amplicon_fh.write(name + "\tCargo\t" + cargo_amplicon + "\n")

            # define quantification window
            cargo_qw1 = "Cargo:spacer_cut:" + str(spacer_info["cut"]) + "-" + str(spacer_info["cut"] + 1) + ":0"
            cargo_qw2 = "Cargo:AttL1:" + str(spacer_info["cut"] + 1) + "-" + str(spacer_info["cut"] + 30) + ":0"

        # Run Crispresso2
        error_fh = open(os.path.join(output, name + ".job.log"), 'wb')

        ## mkdir folder
        preprocess_output = os.path.join(output, name + "_preprocess_output")
        if not os.path.exists(preprocess_output):
            os.mkdir(preprocess_output)

        ## merge R1 and R2
        if gene_strand == "-":
            subprocess.call(
                "flash %s %s --min-overlap 10 --max-overlap 100 --allow-outies -z -d %s" % (R2, R1, preprocess_output),
                shell=True, stderr=error_fh, stdout=error_fh)
        else:
            subprocess.call(
                "flash %s %s --min-overlap 10 --max-overlap 100 --allow-outies -z -d %s" % (R1, R2, preprocess_output),
                shell=True, stderr=error_fh, stdout=error_fh)

        ## align to plasmid
        if pd.notna(sample["Payload ID"]):
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

        ## cs2, quantification, plot
        if assay == "Control":
            if wt_qw1:
                subprocess.call(
                    "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name WT --guide_seq %s "
                    "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                    "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s" % (
                        unmapped_fastq, wt_amplicon, ng_info["seq"], name, output, ncpu), stderr=error_fh,
                    stdout=error_fh,
                    shell=True)
                subprocess.call("python /home/ubuntu/bin/parse_quantification_windows.py -f %s -o %s -qw %s" % (
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
                "python /home/ubuntu/bin/plotCustomAllelePlot.py -f %s -o %s --plot_center %s --plot_left %s "
                "--plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), 0, 1,
                    len(wt_amplicon) - 1), stderr=error_fh, stdout=error_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/allele2html.py -f %s -r %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "WT"), stderr=error_fh, stdout=error_fh, shell=True)

        elif assay == "AmpSeq":
            subprocess.call(
                "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name WT,Beacon --guide_seq %s "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s "
                "--needleman_wunsch_gap_extend 0" % (
                    unmapped_fastq, wt_amplicon + "," + beacon_amplicon, spacer_info["seq"], name, output, ncpu),
                stderr=error_fh, stdout=error_fh, shell=True)
            if type == "single_atg":
                subprocess.call(
                    "python /home/ubuntu/bin/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s -qw %s" % (
                        os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                        wt_qw1, wt_qw2, beacon_qw1, beacon_qw3), stderr=error_fh, stdout=error_fh, shell=True)
            elif type == "dual_atg":
                subprocess.call(
                    "python /home/ubuntu/bin/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s" % (
                        os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                        wt_qw1, wt_qw2, beacon_qw1), stderr=error_fh, stdout=error_fh, shell=True)
            subprocess.call(
                "python /home/ubuntu/bin/plotCustomAllelePlot.py -f %s -o %s -a WT --plot_center %s --plot_left %s "
                "--plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    spacer_info["cut"] - 1, spacer_info["cut"], len(wt_amplicon) - spacer_info["cut"]), stderr=error_fh,
                stdout=error_fh, shell=True)
            subprocess.call(
                "python /home/ubuntu/bin/plotCustomAllelePlot.py -f %s -o %s -a Beacon --plot_center %s --plot_left "
                "%s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    spacer_info["cut"] - 1, spacer_info["cut"], len(beacon_amplicon) - spacer_info["cut"]),
                stderr=error_fh,
                stdout=error_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "WT", wt_qw1), stderr=error_fh, stdout=error_fh,
                            shell=True)
            subprocess.call("python /home/ubuntu/bin/allele2html.py -f %s -r %s -b %s" % (
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
                    RP1_fastq, wt_amplicon + "," + beacon_amplicon, spacer_info["seq"], name + "_WT_Beacon", output,
                    ncpu),
                stderr=error_fh, stdout=error_fh, shell=True)
            subprocess.call(
                "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name Cargo "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s" % (
                    RP2_fastq, cargo_amplicon, name + "_Cargo", output, ncpu), stderr=error_fh, stdout=error_fh,
                shell=True)

            if type == "single_atg":
                subprocess.call(
                    "python /home/ubuntu/bin/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s -qw %s "
                    "-qw %s" % (
                        os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"),
                        os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"), wt_qw1, wt_qw2, beacon_qw1,
                        beacon_qw2, beacon_qw3), stderr=error_fh, stdout=error_fh, shell=True)
            elif type == "dual_atg":
                subprocess.call(
                    "python /home/ubuntu/bin/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s" % (
                        os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"),
                        os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"), wt_qw1, wt_qw2, beacon_qw1),
                    stderr=error_fh, stdout=error_fh, shell=True)

            subprocess.call("python /home/ubuntu/bin/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s" % (
                os.path.join(output, "CRISPResso_on_" + name + "_Cargo"),
                os.path.join(output, "CRISPResso_on_" + name + "_Cargo"), cargo_qw1, cargo_qw2), stderr=error_fh,
                            stdout=error_fh, shell=True)

            subprocess.call(
                "python /home/ubuntu/bin/plotCustomAllelePlot.py -f %s -o %s -a WT --plot_center %s --plot_left %s "
                "--plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"),
                    os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"), spacer_info["cut"] - 1,
                    spacer_info["cut"], len(wt_amplicon) - spacer_info["cut"]), stderr=error_fh, stdout=error_fh,
                shell=True)
            subprocess.call(
                "python /home/ubuntu/bin/plotCustomAllelePlot.py -f %s -o %s -a Beacon --plot_center %s --plot_left "
                "%s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"),
                    os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"), spacer_info["cut"] - 1,
                    spacer_info["cut"], len(beacon_amplicon) - spacer_info["cut"]), stderr=error_fh, stdout=error_fh,
                shell=True)
            subprocess.call(
                "python /home/ubuntu/bin/plotCustomAllelePlot.py -f %s -o %s -a Cargo --plot_center %s --plot_left %s "
                "--plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name + "_Cargo"),
                    os.path.join(output, "CRISPResso_on_" + name + "_Cargo"), spacer_info["cut"] - 1,
                    spacer_info["cut"],
                    len(cargo_amplicon) - spacer_info["cut"]), stderr=error_fh, stdout=error_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"), "WT", wt_qw1), stderr=error_fh,
                            stdout=error_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"), "Beacon", beacon_qw1), stderr=error_fh,
                            stdout=error_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name + "_Cargo"), "Cargo", cargo_qw2), stderr=error_fh,
                            stdout=error_fh, shell=True)
        error_fh.close()
    amplicon_fh.close()


if __name__ == "__main__":
    main()
