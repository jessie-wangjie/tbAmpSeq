"""
Summarize the INDELs in the quantification windows
Based on CRISPResso2
"""

import argparse
import os
import zipfile

import pandas as pd
from CRISPResso2 import CRISPRessoCOREResources
from CRISPResso2 import CRISPRessoShared


def arrstr_to_arr(val):
    return [int(x) for x in val[1:-1].split(",")]


def get_row_around_cut_assymetrical(row, start, end):
    return row['Aligned_Sequence'][start:end], row['Reference_Sequence'][start:end], row['Read_Status'] == 'UNMODIFIED', \
           row['n_deleted'], row['n_inserted'], row['n_mutated'], row['#Reads'], row['%Reads']


def main():
    pd.set_option('display.max_rows', None)

    parser = argparse.ArgumentParser(description="Summarize the INDELs in the quantification windows")
    parser.add_argument("-f", "--CRISPResso2_folder", type=str,
                        help="CRISPResso output folder containing finished analysis", required=True)
    parser.add_argument("-o", "--output_folder", type=str, help="Output folder", required=True)
    parser.add_argument("-qw", "--quantification_windows", type=str,
                        help="Amplicon:Window_name:Window_region:flanking_bp. Bp positions in the amplicon sequence specifying the quantification window, 1-index",
                        action="append")

    args = parser.parse_args()
    cs2_info = CRISPRessoShared.load_crispresso_info(args.CRISPResso2_folder)

    if not cs2_info["running_info"]["args"].write_detailed_allele_table:
        raise Exception('CRISPResso run must be run with the parameter --write_detailed_allele_table')

    z = zipfile.ZipFile(
        os.path.join(args.CRISPResso2_folder, cs2_info["running_info"]["allele_frequency_table_zip_filename"]))
    zf = z.open(cs2_info["running_info"]["allele_frequency_table_filename"])
    df_alleles = pd.read_csv(zf, sep="\t")
    df_alleles["ref_positions"] = df_alleles["ref_positions"].apply(arrstr_to_arr)

    # generate the stats JSON for the result schema
    b_json = {"sample": cs2_info["running_info"]["args"].name}
    b_json["Merged R1R2 Read Num"] = cs2_info["running_info"]["alignment_stats"]["N_TOT_READS"]
    b_json["WT Aligned Read Num"] = cs2_info["results"]["alignment_stats"]["counts_total"]["WT"]
    if len(cs2_info["results"]["ref_names"])>1:
        b_json["Beacon Aligned Read Num"] = cs2_info["results"]["alignment_stats"]["counts_total"][cs2_info["results"]["ref_names"][1]]
        b_json["Aligned Percentage"] = (b_json["WT Aligned Read Num"] + b_json["Beacon Aligned Read Num"]) / b_json[
            "Merged R1R2 Read Num"]
        b_json["WT Aligned Percentage"] = b_json["WT Aligned Read Num"] / (
                b_json["WT Aligned Read Num"] + b_json["Beacon Aligned Read Num"])
        b_json["Beacon Placement Percentage"] = b_json["Beacon Aligned Read Num"] / (
                b_json["WT Aligned Read Num"] + b_json["Beacon Aligned Read Num"])

    qw_stats = []
    for window in args.quantification_windows:
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

        if (ref_name == "Beacon" and qw_name == "beacon_whole") or (ref_name == "PE" and qw_name == "RT_whole"):
            if "indels" in stats:
                b_json["Beacon Indel Read Num"] = stats["indels"]
            else:
                b_json["Beacon Indel Read Num"] = 0
            b_json["Beacon Indel Percentage"] = b_json["Beacon Indel Read Num"] / b_json["Beacon Aligned Read Num"]
            if "substitution" in stats:
                b_json["Beacon Sub Read Num"] = stats["substitution"]
            else:
                b_json["Beacon Sub Read Num"] = 0
            b_json["Beacon Sub Percentage"] = b_json["Beacon Sub Read Num"] / b_json["Beacon Aligned Read Num"]
            gg = df.groupby(["classification", "n_indel"]).sum()
            print(gg)
            print(gg.loc["modified", "#Reads"].cumsum())
            pd.DataFrame(gg.loc["modified", "#Reads"]).to_csv(
                args.output_folder + "/CRISPResso_quantification_of_editing_frequency.beacon.txt",sep="\t")

    pd.DataFrame(qw_stats).to_csv(args.output_folder + "/CRISPResso_quantification_of_editing_frequency.detailed.txt",
                                  sep="\t", header=True, index=False, na_rep=0)

    pd.Series(b_json).to_json(args.output_folder + "/CRISPResso_stats.json")


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
    n_indel = payload["insertion_n"] + payload["deletion_n"]
    if set(include_idx).issubset(payload["deletion_positions"]):
        whole_window_deletion = row["#Reads"]
    else:
        whole_window_deletion = 0
    return {"#Reads": row["#Reads"], "classification": classification, "indels": indels, "insertion": insertion,
            "deletion": deletion, "substitution": substitution, "whole_window_deletion": whole_window_deletion,
            "n_indel": n_indel}


if __name__ == "__main__":
    main()
