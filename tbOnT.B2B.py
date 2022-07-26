#!/usr/bin/env python
"""
On-Target analysis for Amp-seq and 3Primer-Seq
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
    parser.add_argument("-m", help='TB id')
    parser.add_argument("-i", help='Path to fastq')
    parser.add_argument("-p", help='Number of CPUs to use', default=4)
    parser.add_argument("-o", help='Output folder')

    args = parser.parse_args()
    fastq = args.i
    tbid = args.m
    ncpu = int(args.p)
    output = args.o

    # Read from metasheet
    cur.execute(
        "select miseq_sample_name, re1.file_registry_id, re2.file_registry_id from ampseq_sample_metasheet$raw "
        "join registry_entity as re1 on re1.id = aaan_id "
        "join registry_entity as re2 on re2.id = pp_id "
        "where genomics_ampseq_project_queue = %s", [tbid])

    for record in cur:
        name, aaan_id, pp_id = record
        print(record)

        # Get primer information
        cur.execute(
            "select p1.chromosome, p1.start, p2.end, p1.genome_build, target_gene.direction_of_transcription from primer_pair "
            "join primer as p1 on p1.id = primer_pair.forward_primer "
            "join primer as p2 on p2.id = primer_pair.reverse_primer "
            "join target_gene on target_gene.id = p1.gene_or_target_name "
            "where primer_pair.file_registry_id$ = %s", [pp_id])
        target_chr, wt_start, wt_end, genome_build, target_strand = cur.fetchone()

        # reference genome
        genome_fa = "/home/ubuntu/annotation/2bit/" + genome_build + ".2bit"

        # get r1 and r2 fastq
        print(name)
        r1 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R1_*")[0]
        r2 = glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R2_*")[0]

        # WT amplicon
        wt_amplicon = get_seq(genome_fa, target_chr, wt_start, wt_end, target_strand)
        wt_qw1 = ""

        # Beacon amplicon
        # locate atgRNA, ngRNA in the amplicon
        if aaan_id.startswith("AN"):
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
                        "where atg_ng.file_registry_id$ = %s", [aaan_id])
            sp_seq, beacon_seq, ng_seq, rt_coord, atg_seq = cur.fetchone()
            sp1_info = get_cut_site(wt_amplicon, sp_seq)
            ng_info = get_cut_site(wt_amplicon, ng_seq)
            coord = re.match(r"\[(.*), (.*)\]", rt_coord)
            rt_info = get_cut_site(wt_amplicon, atg_seq[int(coord.group(1)) - 1:int(coord.group(2))])
            beacon = get_beacon_seq(beacon_seq, sp1_info["strand"])

            # beacon seq
            beacon_amplicon = wt_amplicon[0:sp1_info["cut"]] + beacon + wt_amplicon[sp1_info["cut"]:]

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
        elif aaan_id.startswith("AA"):
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
                        "where atg_atg.file_registry_id$ = %s", [aaan_id])
            sp1_seq, sp2_seq, beacon1_seq, beacon2_seq = cur.fetchone()

            sp1_info = get_cut_site(wt_amplicon, sp1_seq)
            sp2_info = get_cut_site(wt_amplicon, sp2_seq)

            beacon = get_beacon_seq(beacon1_seq, sp1_info["strand"], beacon2_seq, sp2_info["strand"])

            # beacon seq
            beacon_amplicon = wt_amplicon[0:sp1_info["cut"]] + beacon + wt_amplicon[sp2_info["cut"]:]

            # define quantification window
            # WT amplicon, spacer cutting 2bp
            wt_qw1 = "WT:spacer1_cut:" + str(sp1_info["cut"]) + "-" + str(sp1_info["cut"] + 1) + ":0"
            # WT amplicon, spacer2 cutting 2bp
            wt_qw2 = "WT:spacer2_cut:" + str(sp2_info["cut"]) + "-" + str(sp2_info["cut"] + 1) + ":0"
            # Beacon amplicon, whole beacon insertion, w/ flank 10bp
            beacon_qw1 = "Beacon:beacon_whole:" + str(sp1_info["cut"] + 1) + "-" + str(
                sp1_info["cut"] + len(beacon)) + ":10"

        # Run Crispresso2
        error_fh = open(os.path.join(output, name + ".job.log"), 'wb')

        # mkdir folder
        preprocess_output = os.path.join(output, name + "_preprocess_output")
        if not os.path.exists(preprocess_output):
            os.mkdir(preprocess_output)

        # merge r1 and r2
        if target_strand == "antisense":
            subprocess.call(
                "flash %s %s --min-overlap 10 --max-overlap 100 --allow-outies -z -d %s" % (r2, r1, preprocess_output),
                shell=True, stderr=error_fh, stdout=error_fh)
        else:
            subprocess.call(
                "flash %s %s --min-overlap 10 --max-overlap 100 --allow-outies -z -d %s" % (r1, r2, preprocess_output),
                shell=True, stderr=error_fh, stdout=error_fh)

        unmapped_fastq = os.path.join(preprocess_output, "out.extendedFrags.fastq.gz")

        # cs2, quantification
        if aaan_id.startswith("AN"):
            subprocess.call(
                "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name WT,Beacon --guide_seq %s "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s " % (
                    unmapped_fastq, wt_amplicon + "," + beacon_amplicon, sp1_info["seq"] + "," + ng_info["seq"], name,
                    output, ncpu), stderr=error_fh, stdout=error_fh, shell=True)
            subprocess.call(
                "python utils/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s -qw %s -qw %s" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                    wt_qw1, wt_qw2, beacon_qw1, beacon_qw2, beacon_qw3), stderr=error_fh, stdout=error_fh, shell=True)

        elif aaan_id.startswith("AA"):
            subprocess.call(
                "CRISPResso --fastq_r1 %s --amplicon_seq %s --amplicon_name WT,Beacon --guide_seq %s "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --name %s --output_folder %s "
                "--write_detailed_allele_table --place_report_in_output_folder --n_processes %s " % (
                    unmapped_fastq, wt_amplicon + "," + beacon_amplicon, sp1_info["seq"] + "," + sp2_info["seq"], name,
                    output, ncpu), stderr=error_fh, stdout=error_fh, shell=True)
            subprocess.call("python utils/parse_quantification_windows.py -f %s -o %s -qw %s -qw %s -qw %s" % (
                os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), wt_qw1,
                wt_qw2, beacon_qw1), stderr=error_fh, stdout=error_fh, shell=True)

        # plot
        subprocess.call(
            "python /home/ubuntu/bin/tbOnT/plotCustomAllelePlot.py -f %s -o %s -a WT --plot_center %s --plot_left %s --plot_right %s "
            "--min_freq 0.01 --plot_cut_point" % (
                os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                sp1_info["cut"] - 1, sp1_info["cut"], len(wt_amplicon) - sp1_info["cut"]), stderr=error_fh,
                stdout=error_fh, shell=True)
        subprocess.call(
            "python /home/ubuntu/bin/tbOnT/plotCustomAllelePlot.py -f %s -o %s -a Beacon --plot_center %s --plot_left %s --plot_right "
            "%s --min_freq 0.01 --plot_cut_point" % (
                os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name),
                sp1_info["cut"] - 1, sp1_info["cut"], len(beacon_amplicon) - sp1_info["cut"]), stderr=error_fh,
                stdout=error_fh, shell=True)

        subprocess.call("python /home/ubuntu/bin/tbOnT/allele2html.py -f %s -r %s -b %s" % (
            os.path.join(output, "CRISPResso_on_" + name), "WT", wt_qw1), stderr=error_fh, stdout=error_fh, shell=True)
        subprocess.call("python /home/ubuntu/bin/tbOnT/allele2html.py -f %s -r %s -b %s" % (
            os.path.join(output, "CRISPResso_on_" + name), "Beacon", beacon_qw1), stderr=error_fh, stdout=error_fh,
                        shell=True)

        error_fh.close()


if __name__ == "__main__":
    main()
