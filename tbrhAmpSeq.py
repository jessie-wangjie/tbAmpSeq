#!/usr/bin/env python
"""
On-Target analysis for rhAmp-seq
Meta information from Benchling
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
    ngs_id = re.sub(".*([B|C]TB\d+).*", "\\1", tbid)
    cur.execute("select id, name, email, eln_id from ngs_tracking where file_registry_id$ = %s "
                "union "
                "select id, project_owner, email, sequencing_run_id from custom_tracking where file_registry_id$ = %s", [ngs_id, ngs_id])
    ngs_stats["ngs_tracking"], ngs_stats["experimenter"], ngs_stats["email"], ngs_stats["project_name"] = cur.fetchone()

    # Query sample metasheet information for BTB
    cur.execute("select miseq_sample_name, aa1.file_registry_id, pp1.file_registry_id, aa2.file_registry_id, pp2.file_registry_id, "
                "aa3.file_registry_id, pp3.file_registry_id, aa4.file_registry_id, pp4.file_registry_id, "
                "aa5.file_registry_id, pp5.file_registry_id, aa6.file_registry_id, pp6.file_registry_id,"
                "sample_name, plate, well_position from ampseq_sample_metasheet_v2$raw "
                "left join registry_entity as aa1 on aa1.id = aaanpnsg_id "
                "left join registry_entity as pp1 on pp1.id = pp_id "
                "left join registry_entity as aa2 on aa2.id = aa_id_2 "
                "left join registry_entity as pp2 on pp2.id = pp_id_2 "
                "left join registry_entity as aa3 on aa3.id = aa_id_3 "
                "left join registry_entity as pp3 on pp3.id = pp_id_3 "
                "left join registry_entity as aa4 on aa4.id = payload_id "
                "left join registry_entity as pp4 on pp4.id = pgi_pp_id "
                "left join registry_entity as aa5 on aa5.id = payload_id_2 "
                "left join registry_entity as pp5 on pp5.id = pgi_pp_id_2 "
                "left join registry_entity as aa6 on aa6.id = payload_id_3 "
                "left join registry_entity as pp6 on pp6.id = pgi_pp_id_3 "
                "where genomics_ampseq_project_queue = %s", [tbid])

    if cur.rowcount == 0:
        project_fh.write("This project doesn't exist in the Benchling!\n")

    for record in cur.fetchall():
        meta_stats = {}
        sample_name, aaan_id, pp_id, aa2, pp2, aa3, pp3, payload_id, pgi_pp_id, payload_id_2, pgi_pp_id_2, payload_id_3, pgi_pp_id_3, \
        meta_stats["samplename"], meta_stats["plate"], meta_stats["well"] = record
        meta_stats["miseq_sample_name"] = sample_name
        meta_stats["genomics_ampseq_project_queue"] = tbid

        # select specific sample
        if sample and sample_name != sample:
            continue

        # skip if no sample name
        if not sample_name:
            continue

        # skip if no fastq
        if len(glob.glob(os.path.abspath(fastq) + "/" + sample_name + "_L001_*/*_R1_*")) == 0:
            meta_stats.update(aaanid=aaan_id, ppid=pp_id)
            os.makedirs(os.path.join(output, "CRISPResso_on_" + sample_name), exist_ok=True)
            pd.Series(meta_stats).to_json(os.path.join(output, "CRISPResso_on_" + sample_name, "CRISPResso_quilt_stats.json"))
            continue

        # run all AAs
        aa_list = [[aaan_id, pp_id, payload_id, pgi_pp_id], [aa2, pp2, payload_id_2, pgi_pp_id_2], [aa3, pp3, payload_id_3, pgi_pp_id_3]]
        for aa, pp, payload, pgi_pp in aa_list:
            cs2_stats = {}
            print([sample_name, aa, pp, payload, pgi_pp])

            # Get primer information
            cur.execute(
                "select p1.chromosome, p1.start, p1.end, p2.start, p2.end, p1.genome_build, target_gene.direction_of_transcription from primer_pair "
                "join primer as p1 on p1.id = primer_pair.forward_primer "
                "join primer as p2 on p2.id = primer_pair.reverse_primer "
                "join target_gene on target_gene.id = p1.gene_or_target_name "
                "where primer_pair.file_registry_id$ = %s", [pp])
            chr, p1_start, p1_end, p2_start, p2_end, genome_build, target_strand = cur.fetchone()
            wt_start = min(p1_start, p2_start)
            wt_end = max(p1_end, p2_end)

            # reference genome
            genome_build = re.sub(".*/", "", genome_build)
            genome_fa = "/home/ubuntu/annotation/2bit/" + genome_build + ".2bit"

            # get r1 and r2 fastq
            if (target_strand == "antisense" or target_strand == "-") and (p1_start < p2_start):
                r1 = glob.glob(os.path.abspath(fastq) + "/" + sample_name + "_L001_*/*_R2_*")[0]
                r2 = glob.glob(os.path.abspath(fastq) + "/" + sample_name + "_L001_*/*_R1_*")[0]
            else:
                r1 = glob.glob(os.path.abspath(fastq) + "/" + sample_name + "_L001_*/*_R1_*")[0]
                r2 = glob.glob(os.path.abspath(fastq) + "/" + sample_name + "_L001_*/*_R2_*")[0]

            # sample job log
            name = sample_name + "." + aa
            job_fh = open(os.path.join(output, name + ".job.log"), 'wb')

            # WT amplicon
            wt_amplicon = get_seq(genome_fa, chr, wt_start, wt_end, target_strand)

            read_length = CRISPRessoCORE.get_avg_read_length_fastq(r1)
            if len(wt_amplicon) > read_length * 2 - 4:
                cs2 = "--force_merge_pairs "
            elif len(wt_amplicon) >= read_length * 2 - 9 and len(wt_amplicon) <= read_length * 2 - 4:
                cs2 = "--stringent_flash_merging "
            amplicon_fh.write(name + "\tWT\t" + wt_amplicon + "\n")

            # Beacon amplicon
            # Get spacers information
            cur.execute("select sp1.bases, sp2.bases, beacon1.bases, beacon2.bases, beacon.bases, "
                        "att.left_half, att.central_dinucleotides, att.right_half, spp.file_registry_id from atg_atg "
                        "join modified_rna as m1 on m1.id = atg_atg.atg1 "
                        "join modified_rna as m2 on m2.id = atg_atg.atg2 "
                        "join atgrna as a1 on a1.id = m1.rna "
                        "join atgrna as a2 on a2.id = m2.rna "
                        "join dna_oligo as sp1 on sp1.id = a1.spacer "
                        "join dna_oligo as sp2 on sp2.id = a2.spacer "
                        "join dna_sequence as beacon1 on beacon1.id = a1.beacon "
                        "join dna_sequence as beacon2 on beacon2.id = a2.beacon "
                        "join attachment_sequence as att on att.id = atg_atg.expected_beacon "
                        "join dna_sequence as beacon on beacon.id = atg_atg.expected_beacon "
                        "join registry_entity as spp on spp.id = atg_atg.spacer_pair "
                        "where atg_atg.file_registry_id$ = %s", [aa])
            sp1_seq, sp2_seq, beacon1_seq, beacon2_seq, beacon_seq, beacon_left, beacon_di, beacon_right, cs2_stats["spp_id"] = cur.fetchone()

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
            beacon_qw2 = "Beacon:beacon_fwd:" + str(sp1_info["cut"] + 1) + "-" + str(sp1_info["cut"] + len(beacon1_seq)) + ":10"
            beacon_qw3 = "Beacon:beacon_rev:" + str(sp1_info["cut"] + len(beacon) - len(beacon2_seq) + 1) + "-" + str(
                sp1_info["cut"] + len(beacon)) + ":10"

            # run cargo amplicon
            if payload and pgi_pp:
                # Get primer information
                cur.execute("select p1.id, p1.genome, p2.id, p2.genome from primer_pair "
                            "join primer as p1 on p1.id = primer_pair.forward_primer "
                            "join primer as p2 on p2.id = primer_pair.reverse_primer "
                            "where primer_pair.file_registry_id$ = %s", [pgi_pp])
                p1_id, p1_genome, p2_id, p2_genome = cur.fetchone()

                for id, type in [[p1_id, p1_genome], [p2_id, p2_genome]]:
                    if type != "plasmid":
                        cur.execute("select primer.start, primer.end from primer where primer.id = %s", [id])
                        g_start, g_end = cur.fetchone()
                    else:
                        cur.execute("select seq.bases, primer.start, primer.end, dna_sequence.bases, att.left_half, att.central_dinucleotides, att.right_half from primer "
                                    "join plasmid on plasmid.id = primer.gene_or_target_name "
                                    "join dna_sequence as seq on seq.id = primer.gene_or_target_name "
                                    "join attachment_sequence as att on att.id = plasmid.attachment_sequence "
                                    "join dna_sequence on dna_sequence.id = plasmid.attachment_sequence "
                                    "where primer.id = %s", [id])
                        payload_seq, p_start, p_end, attp_seq, attp_left, attp_di, attp_right = cur.fetchone()

                #
                coord = re.match(r"\[(.*),(.*)\]", beacon_left)
                attb_di_info = get_cut_site(beacon_amplicon, beacon_seq[int(coord.group(1)) - 1:int(coord.group(2)) + len(beacon_di)], 0)
                coord = re.match(r"\[(.*),(.*)\]", attp_left)
                attp_di_info = get_cut_site(payload_seq, attp_seq[int(coord.group(1)) - 1:int(coord.group(2)) + len(attp_di)], 0)

                if (wt_start == g_start and (target_strand == "+" or target_strand == "sense")) or (wt_start == g_end and (target_strand == "-" or target_strand == "antisense")):
                    cargo_amplicon = beacon_amplicon[0:attb_di_info["cut"]] + payload_seq[attp_di_info["cut"]:p_end]
                    coord = re.match(r"\[(.*),(.*)\]", attp_right)
                    cargo_qw1 = "Cargo:attL2:" + str(attb_di_info["cut"] + 1) + "-" + str(attb_di_info["cut"] + int(coord.group(2)) - int(coord.group(1)) + 1) + ":0"
                elif (wt_end == g_end and (target_strand == "+" or target_strand == "sense")) or (wt_start == g_start and (target_strand == "-" or target_strand == "antisense")):
                    cargo_amplicon = payload_seq[p_start - 1:attp_di_info["cut"]] + beacon_amplicon[attb_di_info["cut"]:]
                    coord = re.match(r"\[(.*),(.*)\]", attp_left)
                    cargo_qw1 = "Cargo:attR1:" + str(attp_di_info["cut"] - int(coord.group(2)) + int(coord.group(1)) - len(attp_di) - 1) + "-" + str(
                        attp_di_info["cut"] - 1) + ":0"

                subprocess.call(
                    "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT,Beacon,Cargo --guide_seq %s --name %s --output_folder %s "
                    "--min_frequency_alleles_around_cut_to_plot 0.05 --write_detailed_allele_table --needleman_wunsch_gap_extend 0 --default_min_aln_score 80 "
                    "--trim_sequences --trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:10:0:true "
                    "--place_report_in_output_folder --n_processes %s --bam_output --suppress_report %s " % (
                        r1, r2, wt_amplicon + "," + beacon_amplicon + "," + cargo_amplicon, sp1_info["seq"] + "," + sp2_info["seq"], name, output,
                        ncpu, cs2),
                    stderr=job_fh, stdout=job_fh, shell=True)

                cs2_stats.update(rhampseq_window_quantification(os.path.join(output, "CRISPResso_on_" + name),
                                                       [wt_qw1, wt_qw2, beacon_qw1, beacon_qw2, beacon_qw3, cargo_qw1]))

                subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s" % (
                    os.path.join(output, "CRISPResso_on_" + name), "Cargo", cargo_qw1), stderr=job_fh, stdout=job_fh, shell=True)

            else:
                subprocess.call(
                    "CRISPResso --fastq_r1 %s --fastq_r2 %s --amplicon_seq %s --amplicon_name WT,Beacon --guide_seq %s --name %s --output_folder %s "
                    "--min_frequency_alleles_around_cut_to_plot 0.05 --write_detailed_allele_table --needleman_wunsch_gap_extend 0 "
                    "--trim_sequences  --trimmomatic_options_string ILLUMINACLIP:/home/ubuntu/annotation/fasta/TruSeq_CD.fa:0:90:10:0:true "
                    "--place_report_in_output_folder --n_processes %s --bam_output --suppress_report %s " % (
                        r1, r2, wt_amplicon + "," + beacon_amplicon, sp1_info["seq"] + "," + sp2_info["seq"], name, output, ncpu, cs2),
                        stderr=job_fh, stdout=job_fh, shell=True)

                cs2_stats.update(rhampseq_window_quantification(os.path.join(output, "CRISPResso_on_" + name), [wt_qw1, wt_qw2, beacon_qw1, beacon_qw2, beacon_qw3]))

            pd.concat([pd.Series(cs2_stats), pd.Series(meta_stats), ngs_stats]).to_json(os.path.join(output, "CRISPResso_on_" + name, "CRISPResso_benchling_stats.json"))
            cs2_stats.update(aaanid=aa, ppid=pp)

            if os.path.exists(os.path.join(output, "CRISPResso_on_" + name, "out.extendedFrags.fastq.gz")):
                cs2_stats["total_read_num"] = CRISPRessoCORE.get_n_reads_fastq(r1)
                cs2_stats["merged_r1r2_read_num"] = CRISPRessoCORE.get_n_reads_fastq(os.path.join(output, "CRISPResso_on_" + name, "out.extendedFrags.fastq.gz"))
            pd.concat([pd.Series(meta_stats), pd.Series(cs2_stats)]).to_json(os.path.join(output, "CRISPResso_on_" + name, "CRISPResso_quilt_stats.json"))

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
            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "WT", wt_qw1, wt_qw2), stderr=job_fh, stdout=job_fh, shell=True)
            # subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s -b %s -s %s-%s" % (
            #    os.path.join(output, "CRISPResso_on_" + name), "WT", wt_qw1, wt_qw2, sp1_info["5P"], sp2_info["5P"]), stderr=job_fh,
            #                stdout=job_fh, shell=True)

            subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s" % (
                os.path.join(output, "CRISPResso_on_" + name), "Beacon", beacon_qw1), stderr=job_fh, stdout=job_fh, shell=True)
            # subprocess.call("python /home/ubuntu/bin/tbOnT/utils/allele2html.py -f %s -r %s -b %s -s %s-%s" % (
            #    os.path.join(output, "CRISPResso_on_" + name), "Beacon", beacon_qw1, sp1_info["5P"],
            #    sp1_info["cut"] + len(beacon) + len(sp2_info["seq"]) - 3), stderr=job_fh, stdout=job_fh, shell=True)


            job_fh.close()

    amplicon_fh.close()
    project_fh.close()
    ngs_stats["run start"] = run_start
    ngs_stats["run end"] = str(datetime.datetime.now())
    pd.Series(ngs_stats).to_json(os.path.join(output, "run.json"))


if __name__ == "__main__":
    main()
