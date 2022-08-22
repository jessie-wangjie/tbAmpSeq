#!/usr/bin/env python
"""
On-Target analysis for AmpSeq
Information from Benchling
"""

import argparse
import glob

from utils.base import *
from utils.common_functions import *


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
    info = pd.read_excel(args.m, skiprows=22)
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
            if pd.notna(sample["forward primer"]):
                cur.execute(
                    "select dna_oligo.bases, target_gene.chromosome, target_gene.genome_build, target_gene.direction_of_transcription from primer "
                    "join target_gene on target_gene.id = primer.gene_or_target_name "
                    "join dna_oligo on dna_oligo.id = primer.id "
                    "where primer.file_registry_id$ = %s", [sample["forward primer"]])
                fp_seq, target_chr, genome_build, target_strand = cur.fetchone()

                cur.execute(
                    "select dna_oligo.bases from primer "
                    "join dna_oligo on dna_oligo.id = primer.id "
                    "where primer.file_registry_id$ = %s", [sample["reverse primer"]])
                rp_seq = cur.fetchone()[0]
            else:
                cur.execute(
                    "select target_gene.chromosome, target_gene.genome_build, target_gene.direction_of_transcription from dna_oligo "
                    "join primer on primer.id=dna_oligo.id "
                    "join target_gene on target_gene.id = primer.gene_or_target_name "
                    "where dna_oligo.bases = %s", [sample["Forward Primer Sequence"]])
                target_chr, genome_build, target_strand = cur.fetchone()

            genome_build = re.sub(".*/", "", genome_build)
            reference_index = "/home/ubuntu/annotation/bwa_index/" + genome_build
            fp_info = align_primer(fp_seq, reference_index, target_chr, "CACTCTTTCCCTACACGACGCTCTTCCGATCT")
            rp_info = align_primer(rp_seq, reference_index, target_chr, "GGAGTTCAGACGTGTGCTCTTCCGATCT")

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

        # reference genome
        genome_fa = "/home/ubuntu/annotation/2bit/" + genome_build + ".2bit"

        # get r1 and r2 fastq
        name = sample["MiSeq samplesheet name (no space, no underscore)"]
        print(name)
        if target_strand == "antisense":
            r1 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R2_*")[0]
            r2 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R1_*")[0]
        else:
            r1 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R1_*")[0]
            r2 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R2_*")[0]

        error_fh = open(os.path.join(output, name + ".job.log"), 'wb')

        # WT amplicon
        wt_amplicon = get_seq(genome_fa, target_chr, wt_start, wt_end, target_strand)
        amplicon_fh.write(name + "\tWT\t" + wt_amplicon + "\n")

        sp1_info = {}

        # atgRNA-ngRNA
        if sample["AA/AN/SG/PN ID"].startswith("AN"):
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
                        "where atg_ng.file_registry_id$ = %s", [sample["AA/AN/SG/PN ID"]])
            sp_seq, beacon_seq, ng_seq, rt_coord, atg_seq = cur.fetchone()
            sp1_info = get_cut_site(wt_amplicon, sp_seq)
            ng_info = get_cut_site(wt_amplicon, ng_seq)

            coord = re.match(r"\[(.*), (.*)\]", rt_coord)
            rt_info = get_cut_site(wt_amplicon, atg_seq[int(coord.group(1)) - 1:int(coord.group(2))])

            beacon = get_beacon_seq(beacon_seq, sp1_info["strand"])

            # beacon seq
            beacon_amplicon = wt_amplicon[0:sp1_info["cut"]] + beacon + wt_amplicon[rt_info["3P"] - 1:]
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

            subprocess.call(
                "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT,Beacon --guide_seq %s "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s "
                "--bam_output --needleman_wunsch_gap_extend 0" % (
                    r1, r2, wt_amplicon + "," + beacon_amplicon, sp1_info["seq"] + "," + ng_info["seq"], name, output,
                    ncpu), stderr=error_fh, stdout=error_fh, shell=True)

            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s -qw %s -qw %s" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    wt_qw1, wt_qw2, beacon_qw1, beacon_qw2, beacon_qw3), stderr=error_fh, stdout=error_fh, shell=True)

        # pegRNA-ngRNA
        elif sample["AA/AN/SG/PN ID"].startswith("PN"):
            # get spacer sequences, beacon sequences, ngRNA sequences
            cur.execute("select sp.bases, ng.bases, pegrna.rt_coordinate, peg.bases from peg_ng "
                        "join modified_rna as m1 on m1.id=peg_ng.modified_pegrna "
                        "join modified_rna as m2 on m2.id=peg_ng.modified_ngrna "
                        "join pegrna on pegrna.id=m1.rna "
                        "join ngrna on ngrna.id=m2.rna "
                        "join dna_oligo as sp on sp.id=pegrna.spacer "
                        "join dna_oligo as ng on ng.id=ngrna.spacer "
                        "join dna_sequence as peg on peg.id=pegrna.id "
                        "where peg_ng.file_registry_id$ = %s", [sample["AA/AN/SG/PN ID"]])
            sp_seq, ng_seq, rt_coord, peg_seq = cur.fetchone()
            sp1_info = get_cut_site(wt_amplicon, sp_seq)
            ng_info = get_cut_site(wt_amplicon, ng_seq)

            coord = re.match(r"\[(.*), (.*)\]", rt_coord)
            if sp1_info["strand"] == "+":
                rt_seq = reverse_complement(peg_seq[int(coord.group(1)) - 1:int(coord.group(2))])
            else:
                rt_seq = peg_seq[int(coord.group(1)) - 1:int(coord.group(2))].upper()

            beacon_amplicon = wt_amplicon[0:sp1_info["cut"]] + rt_seq + wt_amplicon[sp1_info["cut"] + len(rt_seq):]
            amplicon_fh.write(name + "\tPE\t" + beacon_amplicon + "\n")

            # define quantification window
            # WT amplicon, spacer cutting 2bp
            wt_qw1 = "WT:spacer_cut:" + str(sp1_info["cut"]) + "-" + str(sp1_info["cut"] + 1) + ":0"

            # WT amplicon, ngRNA cutting 2bp
            wt_qw2 = "WT:ng_cut:" + str(ng_info["cut"]) + "-" + str(ng_info["cut"] + 1) + ":0"

            # Beacon amplicon, whole beacon insertion, w/ flank 10bp
            beacon_qw1 = "PE:RT_whole:" + str(sp1_info["cut"] + 1) + "-" + str(sp1_info["cut"] + len(rt_seq)) + ":0"

            # Beacon amplicon, RT 5'
            beacon_qw2 = "PE:RT_3P:" + str(sp1_info["cut"]) + "-" + str(sp1_info["cut"] + 1) + ":0"
            beacon_qw3 = "PE:RT_5P:" + str(sp1_info["cut"] + len(rt_seq) - 1) + "-" + str(
                sp1_info["cut"] + len(rt_seq) + 1) + ":0"

            # Beacon amplicon, ngRNA cutting 2bp
            ng_info = get_cut_site(beacon_amplicon, ng_seq)
            beacon_qw4 = "PE:ng_cut:" + str(ng_info["cut"]) + "-" + str(ng_info["cut"] + 1) + ":0"

            subprocess.call(
                "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT,PE --guide_seq %s "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s "
                "--bam_output --needleman_wunsch_gap_extend 0" % (
                    r1, r2, wt_amplicon + "," + beacon_amplicon, sp1_info["seq"] + "," + ng_info["seq"], name, output,
                    ncpu), stderr=error_fh, stdout=error_fh, shell=True)

            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s -qw %s -qw %s -qw %s" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    wt_qw1, wt_qw2, beacon_qw1, beacon_qw2, beacon_qw3, beacon_qw4), stderr=error_fh, stdout=error_fh, shell=True)

        # locate atgRNA, atgRNA in the amplicon
        elif sample["AA/AN/SG/PN ID"].startswith("AA"):
            # Get spacers information
            cur.execute("select sp1.bases, sp2.bases, beacon1.bases, beacon2.bases from atg_atg "
                        "join modified_rna as m1 on m1.id = atg_atg.atg1 "
                        "join modified_rna as m2 on m2.id = atg_atg.atg2 "
                        "join atgrna as a1 on a1.id = m1.rna "
                        "join atgrna as a2 on a2.id = m2.rna "
                        "join dna_oligo as sp1 on sp1.id = a1.spacer "
                        "join dna_oligo as sp2 on sp2.id = a2.spacer "
                        "join dna_sequence as beacon1 on beacon1.id = a1.beacon "
                        "join dna_sequence as beacon2 on beacon2.id = a2.beacon "
                        "where atg_atg.file_registry_id$ = %s", [sample["AA/AN/SG/PN ID"]])
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

            subprocess.call(
                "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT,Beacon --guide_seq %s "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s "
                "--bam_output --needleman_wunsch_gap_extend 0" % (
                    r1, r2, wt_amplicon + "," + beacon_amplicon, sp1_info["seq"] + "," + sp2_info["seq"], name, output,
                    ncpu), stderr=error_fh, stdout=error_fh, shell=True)

            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    wt_qw1, wt_qw2, beacon_qw1), stderr=error_fh, stdout=error_fh, shell=True)

        # only sgRNA for indel analysis
        elif sample["AA/AN/SG/PN ID"].startswith("SG"):
            cur.execute("select dna_oligo.bases from sgrna "
                        "join dna_oligo on dna_oligo.id=sgrna.spacer "
                        "where sgrna.file_registry_id$ = %s", [sample["AA/AN/SG/PN ID"]])
            sg_seq = cur.fetchone()[0]
            sp1_info = get_cut_site(wt_amplicon, sg_seq)
            wt_qw1 = "WT:sg_cut:" + str(sp1_info["cut"]) + "-" + str(sp1_info["cut"] + 1) + ":0"

            subprocess.call(
                "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT --guide_seq %s "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s "
                "--bam_output --needleman_wunsch_gap_extend 0" % (
                    r1, r2, wt_amplicon, sp1_info["seq"], name, output, ncpu), stderr=error_fh, stdout=error_fh, shell=True)

            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/parse_quantification_windows.py -f %s -o %s -qw %s" % (
                os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                wt_qw1), stderr=error_fh, stdout=error_fh, shell=True)

        else:
            subprocess.call(
                "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s "
                "--bam_output --needleman_wunsch_gap_extend 0" % (
                    r1, r2, wt_amplicon, name, output, ncpu), stderr=error_fh, stdout=error_fh, shell=True)

        # Plot
        if sp1_info:
            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a WT --plot_center %s "
                "--plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    sp1_info["cut"] - 1, sp1_info["cut"], len(wt_amplicon) - sp1_info["cut"]), stderr=error_fh,
                    stdout=error_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "WT", wt_qw1), stderr=error_fh, stdout=error_fh, shell=True)

        else:
            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a WT --plot_center %s "
                "--plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), 0, 1,
                    len(wt_amplicon) - 1), stderr=error_fh, stdout=error_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "WT"), stderr=error_fh, stdout=error_fh, shell=True)

        if beacon_amplicon:
            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a Beacon --plot_center %s "
                "--plot_left %s --plot_right %s --min_freq 0.01 --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    sp1_info["cut"] - 1, sp1_info["cut"], len(beacon_amplicon) - sp1_info["cut"]), stderr=error_fh,
                    stdout=error_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "Beacon", beacon_qw1), stderr=error_fh, stdout=error_fh, shell=True)

        error_fh.close()
    amplicon_fh.close()


if __name__ == "__main__":
    main()
