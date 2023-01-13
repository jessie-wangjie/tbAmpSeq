#!/usr/bin/env python
"""
Tome 3P-seq analysis
Input from Excel
"""

import argparse
import glob
import os
import re
import subprocess

import pandas as pd


def align_primer(seq, index, adapter=""):
    seq = seq.upper()
    adapter = adapter.upper()
    if adapter in seq:
        seq = seq.replace(adapter, "")

    fastq = ">seq" + "\n" + seq + "\n+" + "\n" + seq
    print(index)
    bwa_out = subprocess.Popen("echo -e '%s' | bwa fastmap %s -" % (fastq, index), stdout=subprocess.PIPE, shell=True)
    pos = subprocess.check_output(["grep", "EM"], stdin=bwa_out.stdout).split()[4:]
    for p in pos:
        if "alt" not in p.decode() and "chr11" not in p.decode():
            print(p.decode())
            m = re.match(r"(.*):([+|-])(\d+)", p.decode())
            break
    return {"chr": m.group(1), "start": int(m.group(3)), "end": int(m.group(3)) + len(seq) - 1, "strand": m.group(2),
            "seq": seq}


def get_cut_site(seq, guide):
    # return 1-index
    # cut is always the left of the cutting site
    guide = guide.replace("U", "T").upper()
    if guide in seq:
        p5 = seq.find(guide) + 1
        p3 = p5 + len(guide) - 1
        cut = p3 - 3
        guide_strand = "+"
    else:
        p3 = seq.find(reverse_complement(guide)) + 1
        p5 = p3 + len(guide) - 1
        cut = p3 - 1 + 3
        guide_strand = "-"
    return {"5P": p5, "3P": p3, "cut": cut, "strand": guide_strand, "seq": guide}


def get_seq(twobit_file, chromosome, start, end, strand):
    seq = subprocess.check_output(
        "twoBitToFa -seq=%s -start=%s -end=%s %s stdout | grep -v \> | xargs | sed 's/ //g'" % (
        chromosome, start - 1, end, twobit_file), shell=True).decode().rstrip()
    if strand == "-":
        seq = reverse_complement(seq)
    return seq.upper()


def get_beacon_seq(seq1, sp1_strand, seq2="", sp2_strand="", attR2="CTCC"):
    beacon = seq1.upper()
    if sp1_strand == "+":
        beacon = reverse_complement(seq1)
    if seq2 != "":
        beacon2 = seq2.upper()
        if sp2_strand == "+":
            beacon2 = reverse_complement(seq2)
        idx = beacon.find(beacon2[0:6])
        beacon = beacon[0:idx] + beacon2
    idx = beacon.find(attR2)
    attL1 = beacon[:idx]
    attR2 = beacon[idx:]
    return {"seq": beacon, "attL1": attL1, "attR2": attR2}


def reverse_complement(seq):
    nt_complement = dict({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '_': '_', '-': '-', 'U': 'A'})
    return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])


def main():
    # Parse command line options
    parser = argparse.ArgumentParser(description='Process meta data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", help='Path to the meta table')
    parser.add_argument("-i", help='Path to fastq')
    parser.add_argument("-p", help='Number of CPUs to use', default=4)
    parser.add_argument("-o", help='Output folder', default=4)
    parser.add_argument("-s", help='specific sample', default=None)
    parser.add_argument("-cs2", help='CRISPRESSO2 parameters', default="")

    args = parser.parse_args()
    fastq = args.i
    ncpu = int(args.p)
    info = pd.read_excel(args.m, skiprows=22)
    output = args.o
    sname = args.s
    cs2 = args.cs2

    amplicon_fh = open(os.path.join(output, os.path.basename(fastq) + ".amplicon.txt"), 'w')

    illumina = {"P5": "CACTCTTTCCCTACACGACGCTCTTCCGATCT", "P7": "GGAGTTCAGACGTGTGCTCTTCCGATCT"}

    for i in info.index:
        sample = info.loc[i]

        if sample["Sample name"] == "PROTOCOLS":
            break

        if sname and sname != sample["MiSeq samplesheet name (no space, no underscore)"]:
            continue

        if pd.isna(sample["Sample name"]):
            continue


        wt_amplicon = ""
        beacon_amplicon = ""
        cargo_amplicon = ""
        assay = "AmpSeq"
        gene_strand = sample["direction of transcription of target gene"]
        name = sample["MiSeq samplesheet name (no space, no underscore)"]

        if pd.notna(sample["attR2"]):
            attR2 = sample["attR2"]
        else:
            attR2 = "CTCC"

        # reference genome
        reference_index = "/home/ubuntu/annotation/bwa_index/" + re.sub(".* ", "", sample["Genome Build"])
        genome_fa = "/home/ubuntu/annotation/2bit/" + re.sub(".* ", "", sample["Genome Build"]) + ".2bit"
        if "Payload ID" in sample and pd.notna(sample["Payload ID"]):
            cargo_fa = "/home/ubuntu/annotation/2bit/" + sample["Payload ID"] + ".2bit"
            if pd.notna(sample["attL2"]):
                attL2_info = align_primer(sample["attL2"],
                                          "/home/ubuntu/annotation/bwa_index/" + sample["Payload ID"])
            else:
                attL2_info = align_primer("CTCAGTGGTGTACGGTACAAA",
                                          "/home/ubuntu/annotation/bwa_index/" + sample["Payload ID"])

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

        if "Reverse Primer 2 Sequence" in sample and pd.notna(sample["Reverse Primer 2 Sequence"]):
            assay = "3P"
            if "Reverse Primer 2 link to Illumina Adapater" not in sample:
                RP2_adapter = "P7"
            else:
                RP2_adapter = sample["Reverse Primer 2 link to Illumina Adapater"]
            RP2_info = align_primer(sample["Reverse Primer 2 Sequence"],
                                    "/home/ubuntu/annotation/bwa_index/" + sample["Payload ID"], illumina[RP2_adapter])
            print(RP2_info)

        print(RP1_info["start"])
        print(FP_info["end"])
        wt_amplicon = get_seq(genome_fa, FP_info["chr"], FP_info["start"], RP1_info["end"], gene_strand)
        #        wt_amplicon=get_seq(genome_fa,FP_info["chr"],RP1_info["start"],FP_info["end"],gene_strand)
        amplicon_fh.write(name + "\tWT\t" + wt_amplicon + "\n")

        spacer_info = {}
        rt_info = {}
        wt_qw1 = ""
        # Beacon amplicon
        ## locate atgRNA, ngRNA in the amplicon
        #        if "atg-ng pair ID" in sample and pd.notna(sample["atg-ng pair ID"]):
        if sample["atgRNA pairing type"] == "single":
            type = "single_atg"
            spacer_info = get_cut_site(wt_amplicon, sample["atgRNA1 spacer sequence"])
            ng_info = get_cut_site(wt_amplicon, sample["atgRNA2 spacer sequence"])
            if pd.notna(sample["atgRNA1 RT seq"]):
                rt_info = get_cut_site(wt_amplicon, sample["atgRNA1 RT seq"])
            beacon = get_beacon_seq(sample["atgRNA1 beacon seq"], spacer_info["strand"], attR2=attR2)
            print(beacon)
            # beacon seq
            beacon_amplicon = wt_amplicon[0:spacer_info["cut"]] + beacon["seq"] + wt_amplicon[spacer_info["cut"]:]
            amplicon_fh.write(name + "\tBeacon\t" + beacon_amplicon + "\n")

            # define quantification window
            # WT amplicon, spacer cutting 2bp
            wt_qw1 = "WT:spacer_cut:" + str(spacer_info["cut"]) + "-" + str(spacer_info["cut"] + 1) + ":0"
            # WT amplicon, ngRNA cutting 2bp
            wt_qw2 = "WT:ng_cut:" + str(ng_info["cut"]) + "-" + str(ng_info["cut"] + 1) + ":0"
            # Beacon amplicon, whole beacon insertion, w/ flank 10bp
            beacon_qw1 = "Beacon:beacon_whole:" + str(spacer_info["cut"] + 1) + "-" + str(
                spacer_info["cut"] + len(beacon["seq"])) + ":10"
            # Beacon amplicon, RT 5'
            if not rt_info:
                beacon_qw2 = ""
            elif rt_info["strand"] == "+":
                beacon_qw2 = "Beacon:RT_5P:" + str(rt_info["5P"] - 1) + "-" + str(rt_info["5P"]) + ":0"
            else:
                beacon_qw2 = "Beacon:RT_5P:" + str(rt_info["5P"]) + "-" + str(rt_info["5P"] + 1) + ":0"
            # Beacon amplicon, ngRNA cutting 2bp
            ng_info = get_cut_site(beacon_amplicon, sample["atgRNA2 spacer sequence"])
            beacon_qw3 = "Beacon:ng_cut:" + str(ng_info["cut"]) + "-" + str(ng_info["cut"] + 1) + ":0"

        ## locate atgRNA, atgRNA in the amplicon
        #        elif "atg-atg pair ID" in sample and pd.notna(sample["atg-atg pair ID"]):
        elif sample["atgRNA pairing type"] == "dual":
            type = "dual_atg"
            spacer_info = get_cut_site(wt_amplicon, sample["atgRNA1 spacer sequence"])
            spacer2_info = get_cut_site(wt_amplicon, sample["atgRNA2 spacer sequence"])
            beacon = get_beacon_seq(sample["atgRNA1 beacon seq"], spacer_info["strand"], sample["atgRNA2 beacon seq"],
                                    spacer2_info["strand"],attR2=attR2)
            print(beacon)
            # beacon seq
            beacon_amplicon = wt_amplicon[0:spacer_info["cut"]] + beacon["seq"] + wt_amplicon[spacer2_info["cut"]:]
            amplicon_fh.write(name + "\tBeacon\t" + beacon_amplicon + "\n")

            # define quantification window
            # WT amplicon, spacer cutting 2bp
            wt_qw1 = "WT:spacer1_cut:" + str(spacer_info["cut"]) + "-" + str(spacer_info["cut"] + 1) + ":0"
            # WT amplicon, spacer2 cutting 2bp
            wt_qw2 = "WT:spacer2_cut:" + str(spacer2_info["cut"]) + "-" + str(spacer2_info["cut"] + 1) + ":0"
            # Beacon amplicon, whole beacon insertion, w/ flank 10bp
            beacon_qw1 = "Beacon:beacon_whole:" + str(spacer_info["cut"] + 1) + "-" + str(
                spacer_info["cut"] + len(beacon["seq"])) + ":10"

        else:
            assay = "Control"
            if "ngRNA ID" in sample and pd.notna(sample["ngRNA ID"]):
                ng_info = get_cut_site(wt_amplicon, sample["ngRNA spacer sequence"])
                wt_qw1 = "WT:ng_cut:" + str(ng_info["cut"]) + "-" + str(ng_info["cut"] + 1) + ":0"

        # cargo seq
        if assay == "3P":
            if sample["Reverse Primer 2 link to Illumina Adapater"] == "P7":
                cargo_amplicon = wt_amplicon[0:spacer_info["cut"]] + beacon["attL1"] + get_seq(cargo_fa,
                                                                                               RP2_info["chr"],
                                                                                               attL2_info["start"],
                                                                                               RP2_info["end"], "+")
                amplicon_fh.write(name + "\tCargo\t" + cargo_amplicon + "\n")

                # define quantification window
                cargo_qw1 = "Cargo:spacer_cut:" + str(spacer_info["cut"]) + "-" + str(spacer_info["cut"] + 1) + ":0"
                cargo_qw2 = "Cargo:AttL1:" + str(spacer_info["cut"] + 1) + "-" + str(spacer_info["cut"] + len(beacon["attL1"])) + ":0"
            else:
                cargo_amplicon = get_seq(cargo_fa, RP2_info["chr"], RP2_info["start"], attL2_info["start"] - 1, "+") + \
                                 beacon["attR2"] + wt_amplicon[spacer2_info["cut"]:]
                amplicon_fh.write(name + "\tCargo\t" + cargo_amplicon + "\n")

                cargo_qw1 = "Cargo:spacer2_cut:" + str(
                    attL2_info["start"] - RP2_info["start"] + len(beacon["attR2"])) + "-" + str(
                    attL2_info["start"] - RP2_info["start"] + len(beacon["attR2"]) + 1) + ":0"
                cargo_qw2 = "Cargo:AttR2:" + str(attL2_info["start"] - RP2_info["start"] + 1) + "-" + str(
                    attL2_info["start"] - RP2_info["start"] + len(beacon["attR2"])) + ":0"

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
#            CRISPRessoShared.force_merge_pairs(os.path.join(preprocess_output, "out.notCombined_1.fastq.gz"),os.path.join(preprocess_output, "out.notCombined_2.fastq.gz"),os.path.join(preprocess_output, "out.forcemerged_uncombined.fastq.gz"))
#            subprocess.call("zcat %s %s | gzip -c > %s" %(os.path.join(preprocess_output, "out.extendedFrags.fastq.gz"), os.path.join(preprocess_output, "out.forcemerged_uncombined.fastq.gz"), os.path.join(preprocess_output, "out.forcemerged.fastq.gz")), shell=True)
#            subprocess.call("cp %s %s" %(os.path.join(preprocess_output, "out.forcemerged.fastq.gz"), os.path.join(preprocess_output, "out.extendedFrags.fastq.gz")), shell=True)
        else:
            subprocess.call(
                "flash %s %s --min-overlap 10 --max-overlap 100 --allow-outies -z -d %s" % (R1, R2, preprocess_output),
                shell=True, stderr=error_fh, stdout=error_fh)

        ## align to plasmid
        if "Payload ID" in sample and pd.notna(sample["Payload ID"]):
            plasmid_left = {"Pdy0186": 140, "PL249": 140, "PL268": 592, "PL210": 15, "PL753": 150, "pDY181": 3155,
                            "PL1113": 402, "OJ611": 616, "AdVG010": 14065, "AdVG013": 14065, "AdVG012": 14065,
                            "ADEG014": 8, "PL983": 1055, "pDY679": 1055, "PL1051": 1055, "PL995": 1055, "PL1178": 3189,
                            "PL1220": 436, "PL1218": 129}
            plasmid_right = {"pDY181": 3195, "PL1113": 442, "OJ611": 656}
            input_file = os.path.join(preprocess_output, "out.extendedFrags.fastq.gz")
            plasmid_sam = os.path.join(preprocess_output, "out.extendedFrags.plasmid.sam")
            unmapped_id = os.path.join(preprocess_output, "out.extendedFrags.unmapped_plasmid.id")
            unmapped_fastq = os.path.join(preprocess_output, "out.extendedFrags.unmapped_plasmid.fastq")
            subprocess.call("bowtie2 --local -p 8 -x %s -U %s -S %s" % (
            "/home/ubuntu/annotation/bowtie2_index/" + sample["Payload ID"], input_file, plasmid_sam), shell=True,
                            stderr=error_fh, stdout=error_fh)
            if sample["Reverse Primer 2 link to Illumina Adapater"] == "P7":
                subprocess.call("samtools view -F4 %s | awk '{ if($4>=%s) print $1}' > %s" % (
                plasmid_sam, plasmid_left[sample["Payload ID"]], unmapped_id), shell=True, stderr=error_fh,
                                stdout=error_fh)
            else:
                subprocess.call("bedtools bamtobed -i %s | awk '{ if($3<=%s) print $4}' > %s" % (
                plasmid_sam, plasmid_right[sample["Payload ID"]], unmapped_id), shell=True, stderr=error_fh,
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
                    "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name WT --guide_seq %s --min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s --write_detailed_allele_table --place_report_in_output_folder --n_processes %s" % (
                    unmapped_fastq, wt_amplicon, ng_info["seq"], name, output, ncpu), stderr=error_fh, stdout=error_fh,
                    shell=True)
                subprocess.call("python /home/ubuntu/bin/parse_quantification_windows.py -f %s -o %s -qw %s" % (
                os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), wt_qw1),
                                stderr=error_fh, stdout=error_fh, shell=True)
            else:
                subprocess.call(
                    "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name WT --min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s --write_detailed_allele_table --place_report_in_output_folder --n_processes %s" % (
                    unmapped_fastq, wt_amplicon, name, output, ncpu), stderr=error_fh, stdout=error_fh, shell=True)

            subprocess.call(
                "python /home/ubuntu/bin/plotCustomAllelePlot.py -f %s -o %s --plot_center %s --plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), 0, 1,
                len(wt_amplicon) - 1), stderr=error_fh, stdout=error_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/allele2html.py -f %s -r %s" % (
            os.path.join(output, "CRISPResso_on_" + name), "WT"), stderr=error_fh, stdout=error_fh, shell=True)

        elif assay == "AmpSeq":
            subprocess.call(
                "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name WT,Beacon --guide_seq %s --min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s --write_detailed_allele_table --place_report_in_output_folder --n_processes %s --needleman_wunsch_gap_extend 0" % (
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
                "python /home/ubuntu/bin/plotCustomAllelePlot.py -f %s -o %s -a WT --plot_center %s --plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                spacer_info["cut"] - 1, spacer_info["cut"], len(wt_amplicon) - spacer_info["cut"]), stderr=error_fh,
                stdout=error_fh, shell=True)
            subprocess.call(
                "python /home/ubuntu/bin/plotCustomAllelePlot.py -f %s -o %s -a Beacon --plot_center %s --plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                spacer_info["cut"] - 1, spacer_info["cut"], len(beacon_amplicon) - spacer_info["cut"]), stderr=error_fh,
                stdout=error_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/allele2html.py -f %s -r %s -b %s" % (
            os.path.join(output, "CRISPResso_on_" + name), "WT", wt_qw1), stderr=error_fh, stdout=error_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/allele2html.py -f %s -r %s -b %s" % (
            os.path.join(output, "CRISPResso_on_" + name), "Beacon", beacon_qw1), stderr=error_fh, stdout=error_fh,
                            shell=True)

        elif assay == "3P":
            ## sort reads by RP1 and RP2, with max 2 mismatches
            RP1_fastq = os.path.join(preprocess_output, "RP1.fastq.gz")
            RP2_fastq = os.path.join(preprocess_output, "RP2.fastq.gz")
            noRP_fastq = os.path.join(preprocess_output, "noRP.fastq.gz")
            if gene_strand == "-":
                if sample["Reverse Primer 2 link to Illumina Adapater"] == "P5":
                    subprocess.call(
                        "/home/ubuntu/software/miniconda3/bin/cutadapt -m 10 -O 10 -e 2.5 -a %s --action=none -o %s --untrimmed-output - %s | /home/ubuntu/software/miniconda3/bin/cutadapt -m 10 -O 10 -e 2.5 -a %s --action=none --untrimmed-output %s -o %s -" % (
                        reverse_complement(FP_info["seq"]), RP1_fastq, unmapped_fastq,
                        reverse_complement(RP2_info["seq"]),
                        noRP_fastq, RP2_fastq), shell=True, stdout=error_fh, stderr=error_fh)
                else:
                    subprocess.call(
                        "/home/ubuntu/software/miniconda3/bin/cutadapt -m 10 -O 10 -e 2.5 -a %s --action=none -o %s --untrimmed-output - %s | /home/ubuntu/software/miniconda3/bin/cutadapt -m 10 -O 10 -e 2.5 -a %s --action=none --untrimmed-output %s -o %s -" % (
                        reverse_complement(FP_info["seq"]), RP1_fastq, unmapped_fastq,
                        reverse_complement(RP2_info["seq"]), noRP_fastq, RP2_fastq), shell=True, stdout=error_fh,
                        stderr=error_fh)
            else:
                if sample["Reverse Primer 2 link to Illumina Adapater"] == "P7":
                    subprocess.call(
                        "/home/ubuntu/software/miniconda3/bin/cutadapt -m 10 -O 10 -e 2.5 -a %s --action=none -o %s --untrimmed-output - %s | /home/ubuntu/software/miniconda3/bin/cutadapt -m 10 -O 10 -e 2.5 -a %s --action=none --untrimmed-output %s -o %s -" % (
                        reverse_complement(RP1_info["seq"]), RP1_fastq, unmapped_fastq,
                        reverse_complement(RP2_info["seq"]), noRP_fastq, RP2_fastq), shell=True, stdout=error_fh,
                        stderr=error_fh)
                else:
                    subprocess.call(
                        "/home/ubuntu/software/miniconda3/bin/cutadapt -m 10 -O 10 -e 2.5 -g %s --action=none -o %s --untrimmed-output - %s | /home/ubuntu/software/miniconda3/bin/cutadapt -m 10 -O 10 -e 2.5 -g %s --action=none --untrimmed-output %s -o %s -" % (
                        FP_info["seq"], RP1_fastq, unmapped_fastq, RP2_info["seq"], noRP_fastq, RP2_fastq), shell=True,
                        stdout=error_fh, stderr=error_fh)
            #            subprocess.call("ls %s" %(RP1_fastq))
            subprocess.call(
                "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name WT,Beacon --guide_seq %s --min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s --write_detailed_allele_table --place_report_in_output_folder --n_processes %s %s" % (
                RP1_fastq, wt_amplicon + "," + beacon_amplicon, spacer_info["seq"], name + "_WT_Beacon", output, ncpu, cs2),
                stderr=error_fh, stdout=error_fh, shell=True)
            subprocess.call(
                "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name Cargo --min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s --write_detailed_allele_table --place_report_in_output_folder --n_processes %s %s" % (
                unmapped_fastq, cargo_amplicon, name + "_Cargo", output, ncpu, cs2), stderr=error_fh, stdout=error_fh, shell=True)

            if type == "single_atg":
                subprocess.call(
                    "python /home/ubuntu/bin/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s -qw %s -qw %s" % (
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
                "python /home/ubuntu/bin/plotCustomAllelePlot.py -f %s -o %s -a WT --plot_center %s --plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"),
                os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"), spacer_info["cut"] - 1,
                spacer_info["cut"], len(wt_amplicon) - spacer_info["cut"]), stderr=error_fh, stdout=error_fh,
                shell=True)
            subprocess.call(
                "python /home/ubuntu/bin/plotCustomAllelePlot.py -f %s -o %s -a Beacon --plot_center %s --plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"),
                os.path.join(output, "CRISPResso_on_" + name + "_WT_Beacon"), spacer_info["cut"] - 1,
                spacer_info["cut"], len(beacon_amplicon) - spacer_info["cut"]), stderr=error_fh, stdout=error_fh,
                shell=True)
            subprocess.call(
                "python /home/ubuntu/bin/plotCustomAllelePlot.py -f %s -o %s -a Cargo --plot_center %s --plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                os.path.join(output, "CRISPResso_on_" + name + "_Cargo"),
                os.path.join(output, "CRISPResso_on_" + name + "_Cargo"), spacer_info["cut"] - 1, spacer_info["cut"],
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
