#!/usr/bin/env python
"""
On-Target analysis for Amp-seq and 3Primer-Seq
Information from Benchling
"""

import argparse
import glob
import datetime
from utils.base import *
from utils.common_functions import *


def main():
    # Parse command line options
    parser = argparse.ArgumentParser(description='Process meta data', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", help='TB id')
    parser.add_argument("-i", help='Path to fastq')
    parser.add_argument("-s", help='Path to fastq', default=None)
    parser.add_argument("-p", help='Number of CPUs to use', default=4)
    parser.add_argument("-o", help='Output folder', default="./")
    parser.add_argument("-cs2", help='CRISPRESSO2 parameters', default="")

    args = parser.parse_args()
    fastq = args.i
    tbid = args.m
    sample = args.s
    ncpu = int(args.p)
    cs2 = args.cs2
    output = args.o

    # create output folder
    os.makedirs(os.path.join(output, "cs2_alignment_html"), exist_ok=True)

    # amplicon information file
    amplicon_fh = open(os.path.join(output, "amplicon.txt"), 'w')
    # project status file
    project_fh = open(os.path.join(output, "status.txt"), 'w')

    # Get information for benchling NGS tracking entity
    run_start = str(datetime.datetime.now())
    if os.path.exists(os.path.join(output, "run.json")):
        ngs_stats = pd.read_json(os.path.join(output, "run.json"), typ="series")
    else:
        ngs_stats = pd.Series(dtype="object")
    ngs_id = re.sub(".*(BTB\d+).*", "\\1", tbid)
    cur.execute("select id, name, email, eln_id from ngs_tracking where file_registry_id$ = %s", [ngs_id])
    ngs_stats["ngs_tracking"], ngs_stats["experimenter"], ngs_stats["email"], ngs_stats["project_name"] = cur.fetchone()

    # Query sample metasheet information for BTB
    cur.execute("select miseq_sample_name, re1.file_registry_id, aaanpnsg_id, re2.file_registry_id, pp_id, "
                "sample_name, mrna_batch_id, modrna_batch_id, primary_cell_lot_id, lnp_batch_id, plate, well_position "
                "from ampseq_sample_metasheet$raw "
                "left join registry_entity as re1 on re1.id = aaanpnsg_id "
                "left join registry_entity as re2 on re2.id = pp_id "
                "where genomics_ampseq_project_queue = %s "
                "union "
                "select miseq_sample_name, re1.file_registry_id, aaanpnsg_id, re2.file_registry_id, pp_id, "
                "sample_name, mrna_batch_id, modrnabatch_id as modrna_batch_id, primary_cell_lot_id, lnp_batch_id, plate, well_position "
                "from ampseq_sample_metasheet_v2$raw "
                "left join registry_entity as re1 on re1.id = aaanpnsg_id "
                "left join registry_entity as re2 on re2.id = pp_id "
                "where genomics_ampseq_project_queue = %s", [tbid, tbid])

    if cur.rowcount == 0:
        project_fh.write("This project doesn't exist in the Benchling!\n")

    for record in cur.fetchall():
        cs2_stats = {}
        name, aaan_id, cs2_stats["aaanid"], pp_id, cs2_stats["ppid"], cs2_stats["samplename"], cs2_stats["mrna_batch_id"], cs2_stats[
            "modatg_batch_id"], cs2_stats["primary_cell_lot_id"], cs2_stats["lnp_batch_id"], cs2_stats["plate"], cs2_stats["well"] = record
        cs2_stats["miseq_sample_name"] = name
        cs2_stats["genomics_ampseq_project_queue"] = tbid
        print([name, aaan_id, pp_id])

        # skip if no sample name or no fastq
        if not name or len(glob.glob(os.path.abspath(fastq) + "/" + name + "_*/*_R1_*")) == 0:
            continue

        # select specific sample
        if sample and name != sample:
            continue

        # Get primer information
        cur.execute("select target_gene.chromosome, p1.start, p2.end, p1.genome_build, target_gene.direction_of_transcription from primer_pair "
                    "join primer as p1 on p1.id = primer_pair.forward_primer "
                    "join primer as p2 on p2.id = primer_pair.reverse_primer "
                    "join target_gene on target_gene.id = p1.gene_or_target_name "
                    "where primer_pair.file_registry_id$ = %s", [pp_id])
        target_chr, wt_start, wt_end, genome_build, target_strand = cur.fetchone()

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
        wt_amplicon = get_seq(genome_fa, target_chr, wt_start, wt_end, target_strand)
        amplicon_fh.write(name + "\tWT\t" + wt_amplicon + "\n")

        sp1_info = {}

        # Beacon amplicon
        # control sample w/o AA/SG/AN/PN id
        if not aaan_id:
            subprocess.call(
                "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT --name %s --output_folder %s --n_processes %s "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --write_detailed_allele_table --needleman_wunsch_gap_extend 0 "
                "--trim_sequences --trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:10:0:true "
                "--bam_output --suppress_report --place_report_in_output_folder %s " % (r1, r2, wt_amplicon, name, output, ncpu, cs2),
                stderr=job_fh, stdout=job_fh, shell=True)

        # atgRNA-ngRNA
        elif aaan_id.startswith("AN"):
            # get spacer sequences, beacon sequences, ngRNA sequences
            cur.execute("select sp.bases, beacon.bases, ng.bases, atgrna.rt_coordinate, atg.bases, spp.file_registry_id from atg_ng "
                        "join modified_rna as m1 on m1.id = atg_ng.atgrna "
                        "join modified_rna as m2 on m2.id = atg_ng.ngrna "
                        "join atgrna on atgrna.id = m1.rna "
                        "join ngrna on ngrna.id = m2.rna "
                        "join dna_oligo as sp on sp.id = atgrna.spacer "
                        "join dna_sequence as beacon on beacon.id = atgrna.beacon "
                        "join dna_oligo as ng on ng.id = ngrna.spacer "
                        "join dna_sequence as atg on atg.id = atgrna.id "
                        "join registry_entity as spp on spp.id = atg_ng.spacer_pair "
                        "where atg_ng.file_registry_id$ = %s", [aaan_id])
            sp_seq, beacon_seq, ng_seq, rt_coord, atg_seq, cs2_stats["spp_id"] = cur.fetchone()
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
                "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT,Beacon --guide_seq %s --name %s --output_folder %s "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --write_detailed_allele_table --needleman_wunsch_gap_extend 0 "
                "--trim_sequences --trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:10:0:true "
                "--place_report_in_output_folder --n_processes %s --bam_output --suppress_report %s " % (
                    r1, r2, wt_amplicon + "," + beacon_amplicon, sp1_info["seq"] + "," + ng_info["seq"], name, output, ncpu, cs2),
                stderr=job_fh, stdout=job_fh, shell=True)

            cs2_stats.update(
                window_quantification(os.path.join(output, "CRISPResso_on_" + name), [wt_qw1, wt_qw2, beacon_qw1, beacon_qw2, beacon_qw3]))

        # atgRNA-atgRNA
        elif aaan_id.startswith("AA"):
            # Get spacers information
            cur.execute("select sp1.bases, sp2.bases, beacon1.bases, beacon2.bases, spp.file_registry_id from atg_atg "
                        "join modified_rna as m1 on m1.id = atg_atg.atg1 "
                        "join modified_rna as m2 on m2.id = atg_atg.atg2 "
                        "join atgrna as a1 on a1.id = m1.rna "
                        "join atgrna as a2 on a2.id = m2.rna "
                        "join dna_oligo as sp1 on sp1.id = a1.spacer "
                        "join dna_oligo as sp2 on sp2.id = a2.spacer "
                        "join dna_sequence as beacon1 on beacon1.id = a1.beacon "
                        "join dna_sequence as beacon2 on beacon2.id = a2.beacon "
                        "join registry_entity as spp on spp.id = atg_atg.spacer_pair "
                        "where atg_atg.file_registry_id$ = %s", [aaan_id])
            sp1_seq, sp2_seq, beacon1_seq, beacon2_seq, cs2_stats["spp_id"] = cur.fetchone()

            sp1_info = get_cut_site(wt_amplicon, sp1_seq)
            sp2_info = get_cut_site(wt_amplicon, sp2_seq)
            beacon = get_beacon_seq(beacon1_seq, sp1_info["strand"], beacon2_seq, sp2_info["strand"])
            if sp1_info["cut"] > sp2_info["cut"]:
                sp1_info = get_cut_site(wt_amplicon, sp2_seq)
                sp2_info = get_cut_site(wt_amplicon, sp1_seq)
                beacon = get_beacon_seq(beacon2_seq, sp1_info["strand"], beacon1_seq, sp2_info["strand"])

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

        # pegRNA-ngRNA
        elif aaan_id.startswith("PN"):
            # get spacer sequences, beacon sequences, ngRNA sequences
            cur.execute(
                "select sp.bases, ng.bases, pegrna.rt_coordinate, pegrna.pbs_coordinate, scaffold.bases, peg.bases, spp.file_registry_id from peg_ng "
                "join modified_rna as m1 on m1.id=peg_ng.modified_pegrna "
                "join modified_rna as m2 on m2.id=peg_ng.modified_ngrna "
                "join pegrna on pegrna.id=m1.rna "
                "join ngrna on ngrna.id=m2.rna "
                "join dna_oligo as sp on sp.id=pegrna.spacer "
                "join dna_oligo as ng on ng.id=ngrna.spacer "
                "join dna_sequence as peg on peg.id=pegrna.id "
                "join dna_sequence as scaffold on scaffold.id=pegrna.scaffold "
                "join registry_entity as spp on spp.id = peg_ng.spacer_pair "
                "where peg_ng.file_registry_id$ = %s", [aaan_id])
            sp_seq, ng_seq, rt_coord, pbs_coord, scaffold_seq, peg_seq, cs2_stats["spp_id"] = cur.fetchone()
            sp1_info = get_cut_site(wt_amplicon, sp_seq)
            ng_info = get_cut_site(wt_amplicon, ng_seq)

            coord = re.match(r"\[(.*), (.*)\]", rt_coord)
            rt_seq = peg_seq[int(coord.group(1)) - 1:int(coord.group(2))].upper()
            coord = re.match(r"\[(.*), (.*)\]", pbs_coord)
            pbs_seq = peg_seq[int(coord.group(1)) - 1:int(coord.group(2))].upper()

            # define quantification window
            # WT amplicon, spacer cutting 2bp
            wt_qw1 = "WT:spacer_cut:" + str(sp1_info["cut"]) + "-" + str(sp1_info["cut"] + 1) + ":0"

            # WT amplicon, ngRNA cutting 2bp
            wt_qw2 = "WT:ng_cut:" + str(ng_info["cut"]) + "-" + str(ng_info["cut"] + 1) + ":0"

            # Beacon amplicon, whole beacon insertion, w/ flank 10bp
            beacon_qw1 = "Prime-edited:RT_whole:" + str(sp1_info["cut"] + 1) + "-" + str(sp1_info["cut"] + len(rt_seq)) + ":0"

            # Beacon amplicon, RT 5'
            beacon_qw2 = "Prime-edited:RT_3P:" + str(sp1_info["cut"]) + "-" + str(sp1_info["cut"] + 1) + ":0"
            beacon_qw3 = "Prime-edited:RT_5P:" + str(sp1_info["cut"] + len(rt_seq)) + "-" + str(sp1_info["cut"] + len(rt_seq) + 1) + ":0"

            subprocess.call(
                "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT --prime_editing_pegRNA_spacer_seq %s "
                "--prime_editing_pegRNA_extension_seq %s --prime_editing_pegRNA_scaffold_seq %s --prime_editing_nicking_guide_seq %s "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --write_detailed_allele_table --needleman_wunsch_gap_extend 0 "
                "--trim_sequences --trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:10:0:true "
                "--name %s --output_folder %s  --place_report_in_output_folder --n_processes %s --suppress_report %s " % (
                    r1, r2, wt_amplicon, sp1_info["seq"], rt_seq + pbs_seq, scaffold_seq, ng_info["seq"], name, output, ncpu, cs2), stderr=job_fh,
                stdout=job_fh, shell=True)

            # Beacon amplicon, ngRNA cutting 2bp
            beacon_amplicon = read_ref_cs2(os.path.join(output, "CRISPResso_on_" + name), "Prime-edited")
            amplicon_fh.write(name + "\tPrime-edited\t" + beacon_amplicon + "\n")

            ng_info = get_cut_site(beacon_amplicon, ng_seq)
            beacon_qw4 = "Prime-edited:ng_cut:" + str(ng_info["cut"]) + "-" + str(ng_info["cut"] + 1) + ":0"

            cs2_stats.update(window_quantification(os.path.join(output, "CRISPResso_on_" + name),
                                                   [wt_qw1, wt_qw2, beacon_qw1, beacon_qw2, beacon_qw3, beacon_qw4]))

        # sgRNA
        elif aaan_id.startswith("SG"):
            cur.execute("select dna_oligo.bases from sgrna "
                        "join dna_oligo on dna_oligo.id=sgrna.spacer "
                        "where sgrna.file_registry_id$ = %s", [aaan_id])
            sg_seq = cur.fetchone()[0]
            sp1_info = get_cut_site(wt_amplicon, sg_seq)
            wt_qw1 = "WT:sg_cut:" + str(sp1_info["cut"]) + "-" + str(sp1_info["cut"] + 1) + ":0"

            subprocess.call(
                "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT --guide_seq %s --name %s --output_folder %s "
                "--min_frequency_alleles_around_cut_to_plot 0.05 --write_detailed_allele_table --needleman_wunsch_gap_extend 0 "
                "--trim_sequences --trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:10:0:true "
                "--place_report_in_output_folder --n_processes %s --bam_output --suppress_report %s " % (
                    r1, r2, wt_amplicon, sp1_info["seq"], name, output, ncpu, cs2), stderr=job_fh, stdout=job_fh, shell=True)

            cs2_stats.update(window_quantification(os.path.join(output, "CRISPResso_on_" + name), [wt_qw1]))

        pd.concat([pd.Series(cs2_stats), ngs_stats]).to_json(os.path.join(output, "CRISPResso_on_" + name, "CRISPResso_benchling_stats.json"))
        cs2_stats.update(aaanid=aaan_id, ppid=pp_id)
        pd.Series(cs2_stats).to_json(os.path.join(output, "CRISPResso_on_" + name, "CRISPResso_quilt_stats.json"))

        # write sample status
        with open(os.path.join(output, "CRISPResso_on_" + name, "CRISPResso_status.txt"), 'r') as f:
            status = f.readline()
            for line in f:
                if "ERROR" in line:
                    status = line
                    p = subprocess.Popen("samtools view %s | cut -f10 | sort | uniq -c | sort -n | tail -1" % (
                        os.path.join(output, "CRISPResso_on_" + name, "CRISPResso_output.bam")), stdout=subprocess.PIPE, shell=True)
                    status = status.rstrip() + "\t" + p.communicate()[0].decode('utf-8').lstrip()
            project_fh.write(name + "\t" + status)

        # plot
        if sp1_info:
            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a WT --min_freq 0.01 "
                "--plot_center %s --plot_left %s --plot_right %s --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), sp1_info["cut"] - 1,
                    sp1_info["cut"],
                    len(wt_amplicon) - sp1_info["cut"]), stderr=job_fh, stdout=job_fh, shell=True)

            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "WT", wt_qw1), stderr=job_fh, stdout=job_fh, shell=True)
        else:
            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a WT --min_freq 0.01 "
                "--plot_center %s --plot_left %s --plot_right %s --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), 0, 1, len(wt_amplicon) - 1),
                stderr=job_fh, stdout=job_fh, shell=True)

            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s" % (os.path.join(output, "CRISPResso_on_" + name), "WT"),
                            stderr=job_fh, stdout=job_fh, shell=True)

        if aaan_id and aaan_id.startswith("PN"):
            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a Prime-edited --min_freq 0.01 "
                "--plot_center %s --plot_left %s --plot_right %s --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), sp1_info["cut"] - 1,
                    sp1_info["cut"],
                    len(beacon_amplicon) - sp1_info["cut"]), stderr=job_fh, stdout=job_fh, shell=True)
            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a Scaffold-incorporated --min_freq 0.01 "
                "--plot_center %s --plot_left %s --plot_right %s --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), sp1_info["cut"] - 1,
                    sp1_info["cut"],
                    len(beacon_amplicon) - sp1_info["cut"]), stderr=job_fh, stdout=job_fh, shell=True)

            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "Prime-edited", beacon_qw1), stderr=job_fh, stdout=job_fh, shell=True)
            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "Scaffold-incorporated", beacon_qw1), stderr=job_fh, stdout=job_fh, shell=True)

        elif aaan_id and (not aaan_id.startswith("SG")):
            subprocess.call(
                "python /home/ubuntu/bin/tbOnT/utils/plotCustomAllelePlot.py -f %s -o %s -a Beacon --min_freq 0.01 "
                "--plot_center %s --plot_left %s --plot_right %s --plot_cut_point" % (
                    os.path.join(output, "CRISPResso_on_" + name), os.path.join(output, "CRISPResso_on_" + name), sp1_info["cut"] - 1,
                    sp1_info["cut"],
                    len(beacon_amplicon) - sp1_info["cut"]), stderr=job_fh, stdout=job_fh, shell=True)

            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "Beacon", beacon_qw1), stderr=job_fh, stdout=job_fh, shell=True)

        job_fh.close()

    amplicon_fh.close()
    project_fh.close()
    ngs_stats["run start"] = run_start
    ngs_stats["run end"] = str(datetime.datetime.now())
    pd.Series(ngs_stats).to_json(os.path.join(output, "run.json"))


if __name__ == "__main__":
    main()
