#!/usr/bin/env python
"""
Tome 3P-seq analysis
Input from Excel
"""

import argparse
import glob

from utils.common_functions import *


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

        if sname and sname != sample["Sample name"]:
            continue

        if pd.isna(sample["Sample name"]):
            continue

        cargo_amplicon = ""
        gene_strand = sample["direction of transcription of target gene"]
        name = sample["MiSeq samplesheet name (no space, no underscore)"]

        # reference genome
        reference_index = "/home/ubuntu/annotation/bwa_index/" + re.sub(".* ", "", sample["Genome Build"])
        genome_fa = "/home/ubuntu/annotation/2bit/" + re.sub(".* ", "", sample["Genome Build"]) + ".2bit"
        if "Payload ID" in sample and pd.notna(sample["Payload ID"]):
            cargo_fa = "/home/ubuntu/annotation/2bit/" + sample["Payload ID"] + ".2bit"

        # get R1 and R2
        print(name)
        R1 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R1_*")[0]
        R2 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R2_*")[0]

        # locate the FP and RP in the genome
        if "Genomic Primer link to Illumina adapter" not in sample:
            FP_adapter = "P5"
        else:
            FP_adapter = sample["Genomic Primer link to Illumina adapter"]
        FP_info = align_primer(sample["Genomic Primer Sequence"], reference_index, sample["cryptic chromosome"], illumina[FP_adapter])
        print(FP_info)

        if "Cargo Primer link to Illumina adapater" not in sample:
            RP_adapter = "P7"
        else:
            RP_adapter = sample["Cargo Primer link to Illumina adapater"]
        RP_info = align_primer(sample["Cargo Primer Sequence"],
                               "/home/ubuntu/annotation/bwa_index/" + sample["Payload ID"], sample["Payload ID"], illumina[RP_adapter])

        # cargo seq
        cut = int(sample["cryptic coordinate"])
        attL2_info = align_primer("CTCAGTGGTGTACGGTACAAA", "/home/ubuntu/annotation/bwa_index/" + sample["Payload ID"], sample["Payload ID"])
        attR1_info = align_primer("TTTGTCTGGTCAACCACCGCGGT", "/home/ubuntu/annotation/bwa_index/" + sample["Payload ID"], sample["Payload ID"])

        if gene_strand == "-":
            if FP_info["start"] < cut:
                cargo_amplicon = get_seq(cargo_fa, RP_info["chr"], RP_info["start"], attL2_info["start"] - 1,
                                         "+") + get_seq(genome_fa, FP_info["chr"], cut, FP_info["end"], gene_strand)
                cargo_qw1 = "Cargo:AttR1:" + str(
                    attL2_info["start"] - RP_info["start"] - len(attR1_info["seq"])) + "-" + str(
                        attL2_info["start"] - RP_info["start"]) + ":0"
            else:
                cargo_amplicon = get_seq(genome_fa, FP_info["chr"], FP_info["start"], cut - 1, gene_strand) + get_seq(cargo_fa,
                                         RP_info["chr"], attL2_info["start"], RP_info["end"], "+")
                cargo_qw1 = "Cargo:AttL2:" + str(cut - FP_info["start"]) + "-" + str(cut - FP_info["start"] + len(attL2_info["seq"])) + ":0"
        else:
            if FP_info["start"] < cut:
                cargo_amplicon = get_seq(genome_fa, FP_info["chr"], FP_info["start"], cut - 1, gene_strand) + get_seq(cargo_fa,
                                         RP_info["chr"], attL2_info[ "start"], RP_info["end"], "+")
                cargo_qw1 = "Cargo:AttL2:" + str(cut - FP_info["start"]) + "-" + str(cut - FP_info["start"] + len(attL2_info["seq"])) + ":0"
            else:
                cargo_amplicon = get_seq(cargo_fa, RP_info["chr"], RP_info["start"], attL2_info["start"] - 1,
                                         "+") + get_seq(genome_fa, FP_info["chr"], cut, FP_info["end"], gene_strand)
                cargo_qw1 = "Cargo:AttR1:" + str(
                    attL2_info["start"] - RP_info["start"] - len(attR1_info["seq"])) + "-" + str(
                        attL2_info["start"] - RP_info["start"]) + ":0"

        amplicon_fh.write(name + "\tCargo\t" + cargo_amplicon + "\n")

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
        if "Payload ID" in sample and pd.notna(sample["Payload ID"]):
            plasmid_left = {"Pdy0186": 140, "PL249": 140, "PL268": 592, "PL210": 15, "PL753": 150, "pDY181": 3155,
                            "PL1113": 402, "OJ611": 616, "AdVG010": 14065, "AdVG013": 14065, "AdVG012": 14065,
                            "ADEG014": 8, "PL983": 1055, "pDY679": 1055, "PL1051": 1055, "PL995": 1055}
            plasmid_right = {"pDY181": 3195, "PL1113": 442, "OJ611": 656}
            input_file = os.path.join(preprocess_output, "out.extendedFrags.fastq.gz")
            plasmid_sam = os.path.join(preprocess_output, "out.extendedFrags.plasmid.sam")
            unmapped_id = os.path.join(preprocess_output, "out.extendedFrags.unmapped_plasmid.id")
            unmapped_fastq = os.path.join(preprocess_output, "out.extendedFrags.unmapped_plasmid.fastq")
            subprocess.call("bowtie2 --local -p 8 -x %s -U %s -S %s" % (
                "/home/ubuntu/annotation/bowtie2_index/" + sample["Payload ID"], input_file, plasmid_sam), shell=True,
                            stderr=error_fh, stdout=error_fh)
            if sample["Cargo Primer link to Illumina Adapater"] == "P7":
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
        subprocess.call(
            "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name Cargo --min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s --write_detailed_allele_table --place_report_in_output_folder --n_processes %s %s" % (
                unmapped_fastq, cargo_amplicon, name, output, ncpu, cs2), stderr=error_fh, stdout=error_fh, shell=True)

        subprocess.call("python /home/ubuntu/bin/parse_quantification_windows.py -f %s -o %s -qw %s" % (
            os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), cargo_qw1),
                        stderr=error_fh, stdout=error_fh, shell=True)

        subprocess.call(
            "python /home/ubuntu/bin/plotCustomAllelePlot.py -f %s -o %s -a Cargo --plot_center %s --plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                0, 1, len(cargo_amplicon) - 1), stderr=error_fh, stdout=error_fh, shell=True)
        subprocess.call("python /home/ubuntu/bin/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "Cargo", cargo_qw1), stderr=error_fh, stdout=error_fh,
                        shell=True)

        error_fh.close()
    amplicon_fh.close()


if __name__ == "__main__":
    main()
