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


def get_modified_in_quantification_window(row, include_idx):
    start = row['ref_positions'].index(min(include_idx))
    end = row['ref_positions'].index(max(include_idx))
    all_deletion_positions = set(eval(row['all_deletion_positions'])) & include_idx
    all_substitution_positions = set(eval(row['all_substitution_positions'])) & include_idx

    payload = CRISPRessoCOREResources.find_indels_substitutions(row["Aligned_Sequence"], row["Reference_Sequence"], include_idx)
    payload["deletion_n"] = len(all_deletion_positions)
    classification = "unmodified"
    if payload["insertion_n"] + payload["deletion_n"] + payload["substitution_n"] > 0:
        classification = "modified"
    insertion = row["#Reads"] * (payload["insertion_n"] > 0)
    deletion = row["#Reads"] * (payload["deletion_n"] > 0)
    substitution = row["#Reads"] * (payload["substitution_n"] > 0)
    indels = row["#Reads"] * ((payload["insertion_n"] > 0) | (payload["deletion_n"] > 0))
    n_indel = payload["insertion_n"] + payload["deletion_n"]
    n_modified = payload["insertion_n"] + payload["deletion_n"] + payload["substitution_n"]
    if set(include_idx).issubset(payload["deletion_positions"]):
        whole_window_deletion = row["#Reads"]
    else:
        whole_window_deletion = 0
    return {"Aligned_Sequence": row['Aligned_Sequence'][start:end + 1], "Reference_Sequence": row['Reference_Sequence'][start:end + 1],
            "ref_positions": row['ref_positions'][start:end + 1],
            "all_deletion_positions": list(all_deletion_positions), "all_substitution_positions": list(all_substitution_positions),
            "#Reads": row["#Reads"], "classification": classification, "indels": indels, "insertion": insertion, "deletion": deletion,
            "substitution": substitution, "whole_window_deletion": whole_window_deletion, "n_indel": n_indel, "n_modified": n_modified,
            "n_insertion": payload["insertion_n"], "n_deletion": payload["deletion_n"], "n_substitution": payload["substitution_n"]}


def main():
    pd.set_option('display.max_rows', None)

    parser = argparse.ArgumentParser(description="Summarize the INDELs in the quantification windows")
    parser.add_argument("-f", "--CRISPResso2_folder", type=str, help="CRISPResso output folder containing finished analysis", required=True)
    parser.add_argument("-o", "--output_folder", type=str, help="Output folder", required=True)
    parser.add_argument("-qw", "--quantification_windows", type=str, action="append",
                        help="Amplicon:Window_name:Window_region:flanking_bp. Bp positions in the amplicon sequence specifying the quantification window, 1-index")

    args = parser.parse_args()
    # cs2_info = CRISPRessoShared.load_crispresso_info(args.CRISPResso2_folder)

    # if not cs2_info["running_info"]["args"].write_detailed_allele_table:
    #    raise Exception('CRISPResso run must be run with the parameter --write_detailed_allele_table')

    # z = zipfile.ZipFile(os.path.join(args.CRISPResso2_folder, cs2_info["running_info"]["allele_frequency_table_zip_filename"]))
    # zf = z.open(cs2_info["running_info"]["allele_frequency_table_filename"])
    zf = "CRISPResso_on_10/Alleles_frequency_table.zip"
    df_alleles = pd.read_csv(zf, sep="\t")
    df_alleles["ref_positions"] = df_alleles["ref_positions"].apply(arrstr_to_arr)

    for window in args.quantification_windows:
        ref_name, qw_name, qw, flank_bp = window.split(":")
        start, end = qw.split("-")

        df_ref = df_alleles[df_alleles["Reference_Name"] == ref_name]
        if df_ref.empty:
            continue

        df = df_ref.apply(lambda row: get_modified_in_quantification_window(row, set(range(int(start) - 1, int(end)))), axis=1, result_type='expand')
        df = df.groupby(["Aligned_Sequence", "Reference_Sequence", "classification", "n_indel", "n_modified", "n_insertion", "n_deletion",
                         "n_substitution"]).aggregate(
            {"ref_positions": "first", "all_deletion_positions": "first", "all_substitution_positions": "first", "#Reads": "sum", "indels": "sum",
             "insertion": "sum", "deletion": "sum", "substitution": "sum", "whole_window_deletion": "sum"}).reset_index().set_index(
            'Aligned_Sequence').sort_values(["n_modified", "n_indel"])
        df.to_csv(args.output_folder + "/CRISPResso_beacon.alleles_frequency.txt", sep="\t", header=True, na_rep=0)

if __name__ == "__main__":
    main()
