#!/usr/bin/env python
"""
On-Target analysis for Amp-seq and 3Primer-Seq
Information from Benchling

functions need by tbOnT pipeline
"""

import os
import re
import subprocess
import zipfile

import pandas as pd
from CRISPResso2 import CRISPRessoCORE
from CRISPResso2 import CRISPRessoCOREResources
from CRISPResso2 import CRISPRessoShared


def align_primer(seq, index, chromosome, adapter=""):
    seq = seq.upper()
    adapter = adapter.upper()

    if adapter in seq:
        seq = seq.replace(adapter, "")

    fastq = ">seq" + "\n" + seq + "\n+" + "\n" + seq
    bwa_out = subprocess.Popen("echo -e '%s' | bwa fastmap %s -" % (fastq, index), stdout=subprocess.PIPE, shell=True)
    pos = subprocess.check_output(["grep", "EM"], stdin=bwa_out.stdout).split()[4:]
    for p in pos:
        m = re.match(r"(.*):([+|-])(\d+)", p.decode())
        if chromosome == m.group(1):
            return {"chr": m.group(1), "start": int(m.group(3)), "end": int(m.group(3)) + len(seq) - 1,
                    "strand": m.group(2), "seq": seq}


def get_cut_site(seq, guide):
    # return 1-index
    # cut is always the left of the cutting site
    guide = guide.upper().replace("U", "T")
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
    if strand == "-" or strand == "antisense":
        seq = reverse_complement(seq)
    return seq.upper()


def get_beacon_seq(seq1, sp1_strand, seq2="", sp2_strand=""):
    # beacon1 and beacon2 should have at least 5nt overlap.
    beacon = seq1.upper()
    if sp1_strand == "+":
        beacon = reverse_complement(seq1)
    if seq2 != "":
        beacon2 = seq2.upper()
        if sp2_strand == "+":
            beacon2 = reverse_complement(seq2)
        idx = beacon.find(beacon2[0:5])
        beacon = beacon[0:idx] + beacon2
    return beacon


def reverse_complement(seq):
    nt_complement = dict({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '_': '_', '-': '-', 'U': 'A'})
    return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])


def arrstr_to_arr(val):
    return [int(x) for x in val[1:-1].split(",")]


def get_row_around_cut_assymetrical(row, start, end):
    return row['Aligned_Sequence'][start:end], row['Reference_Sequence'][start:end], row['Read_Status'] == 'UNMODIFIED', \
           row['n_deleted'], row['n_inserted'], row['n_mutated'], row['#Reads'], row['%Reads']


def get_modified_in_quantification_window(row, include_idx):
    payload = CRISPRessoCOREResources.find_indels_substitutions(row["Aligned_Sequence"], row["Reference_Sequence"],
                                                                include_idx)
    classification = "unmodified"
    if payload["insertion_n"] + payload["deletion_n"] + payload["substitution_n"] > 0:
        classification = "modified"
    insertion = row["#Reads"] * (payload["insertion_n"] > 0)
    deletion = row["#Reads"] * (payload["deletion_n"] > 0)
    substitution = row["#Reads"] * (payload["substitution_n"] > 0)
    indels = row["#Reads"] * ((payload["insertion_n"] > 0) | (payload["deletion_n"] > 0))
    if set(include_idx).issubset(payload["deletion_positions"]):
        whole_window_deletion = row["#Reads"]
    else:
        whole_window_deletion = 0
    return {"#Reads": row["#Reads"], "classification": classification, "indels": indels, "insertion": insertion,
            "deletion": deletion, "substitution": substitution, "whole_window_deletion": whole_window_deletion}


def window_quantification(cs2_folder, quantification_windows):
    # Amplicon:Window_name:Window_region:flanking_bp. 
    # Bp positions in the amplicon sequence specifying the quantification window, 1-index

    try:
        cs2_info = CRISPRessoShared.load_crispresso_info(cs2_folder)
    except Exception:
        return {}

    if not cs2_info["running_info"]["args"].write_detailed_allele_table:
        raise Exception('CRISPResso run must be run with the parameter --write_detailed_allele_table')

    z = zipfile.ZipFile(os.path.join(cs2_folder, cs2_info["running_info"]["allele_frequency_table_zip_filename"]))
    zf = z.open(cs2_info["running_info"]["allele_frequency_table_filename"])
    df_alleles = pd.read_csv(zf, sep="\t")
    df_alleles["ref_positions"] = df_alleles["ref_positions"].apply(arrstr_to_arr)

    # generate the stats JSON for the result schema
    b_json = {"total_read_num": CRISPRessoCORE.get_n_reads_fastq(cs2_info["running_info"]["args"].fastq_r1),
              "merged_r1r2_read_num": int(cs2_info["running_info"]["alignment_stats"]["N_TOT_READS"]),
              "wt_aligned_read_num": int(cs2_info["results"]["alignment_stats"]["counts_total"]["WT"])}
    if "Beacon" in cs2_info["results"]["alignment_stats"]["counts_total"]:
        b_json["beacon_aligned_read_num"] = int(cs2_info["results"]["alignment_stats"]["counts_total"]["Beacon"])
    elif "PE" in cs2_info["results"]["alignment_stats"]["counts_total"]:
        b_json["beacon_aligned_read_num"] = int(cs2_info["results"]["alignment_stats"]["counts_total"]["PE"])
    else:
        b_json["beacon_aligned_read_num"] = 0

    b_json["aligned_percentage"] = format(100 * (b_json["wt_aligned_read_num"] + b_json["beacon_aligned_read_num"]) / b_json[
        "merged_r1r2_read_num"], ".1f")
    b_json["wt_aligned_percentage"] = format(100 * b_json["wt_aligned_read_num"] / (
        b_json["wt_aligned_read_num"] + b_json["beacon_aligned_read_num"]), ".1f")
    b_json["beacon_placement_percentage"] = format(100 * b_json["beacon_aligned_read_num"] / (
        b_json["wt_aligned_read_num"] + b_json["beacon_aligned_read_num"]), ".1f")

    qw_stats = []
    for window in quantification_windows:
        print(window)
        ref_name, qw_name, qw, flank_bp = window.split(":")
        start, end = qw.split("-")

        stats = {"amplicon": ref_name, "window_name": qw_name, "window_region": qw + ":" + flank_bp}
        df_ref = df_alleles[df_alleles["Reference_Name"] == ref_name]
        if df_ref.empty:
            continue

        df = df_ref.apply(lambda row: get_modified_in_quantification_window(row, set(range(int(start) - 1, int(end)))),
                          axis=1, result_type='expand')
        g = df.groupby("classification").sum()
        for i in g.index:
            stats[i] = g.loc[i]["#Reads"]
            if i == "modified":
                stats.update(g.loc[i][["indels", "insertion", "deletion", "substitution", "whole_window_deletion"]])

        if int(flank_bp):
            include_idx = []
            include_idx.extend(range(int(start) - int(flank_bp) - 1, int(start) - 1))
            include_idx.extend(range(int(end), int(end) + int(flank_bp)))
            df_flank = pd.concat([df, df_ref.apply(
                lambda row: get_modified_in_quantification_window(row, sorted(include_idx)), axis=1,
                result_type='expand').add_suffix("_flank").drop("#Reads_flank", axis=1)], axis=1)
            g = df_flank.groupby(["classification", "classification_flank"]).sum()
            for i, j in g.index:
                stats[i + "_" + j + "_flank"] = g.loc[i, j]["#Reads"]
                if j == "modified":
                    stats.update(g.loc[i, j][g.columns.str.endswith("_flank")].add_prefix(i + "_"))

        qw_stats.append(stats)

        if ref_name == "Beacon" and qw_name == "beacon_whole":
            if "indels" in stats:
                b_json["beacon_indel_read_num"] = int(stats["indels"])
                b_json["beacon_sub_read_num"] = int(stats["substitution"])
            else:
                b_json["beacon_indel_read_num"] = 0
                b_json["beacon_sub_read_num"] = 0
            b_json["beacon_indel_percentage"] = format(100 * b_json["beacon_indel_read_num"] / b_json["beacon_aligned_read_num"], ".1f")
            b_json["beacon_sub_percentage"] = format(100 * b_json["beacon_sub_read_num"] / b_json["beacon_aligned_read_num"], ".1f")
            b_json["perfect_beacon_percent"] = format(
                100 * (b_json["beacon_aligned_read_num"] - b_json["beacon_indel_read_num"]) / (
                        b_json["wt_aligned_read_num"] + b_json["beacon_aligned_read_num"]), ".1f")

    pd.DataFrame(qw_stats).to_csv(cs2_folder + "/CRISPResso_quantification_of_editing_frequency.detailed.txt",
                                  sep="\t", header=True, index=False, na_rep=0)
    pd.Series(b_json).to_json(cs2_folder + "/CRISPResso_stats.json")

    return b_json


def read_ref_cs2(cs2_folder, ref_name):
    try:
        cs2_info = CRISPRessoShared.load_crispresso_info(cs2_folder)
        return cs2_info["results"]["refs"][ref_name]["sequence"]
    except Exception:
        return

