"""
Summarize the scaffold read through
Based on CRISPResso2
"""
import argparse
import os.path
import zipfile
import pandas as pd
from CRISPResso2 import CRISPRessoShared
from base import *
from common_functions import *


def arrstr_to_arr(val):
    return [int(x) for x in val[1:-1].split(",")]


def get_scaffold_search(amplicon, extension_seq, scaffold_seq, scaffold_min_match_length):
    extension_seq = extension_seq.upper().replace('U', 'T')
    len_scaffold_to_use = scaffold_min_match_length

    if extension_seq in amplicon:
        print(extension_seq)
        scaffold_start_loc = amplicon.index(extension_seq) + len(extension_seq)
        print(scaffold_start_loc)
        scaffold_dna_search = scaffold_seq[-len_scaffold_to_use:] + extension_seq
        while scaffold_dna_search in amplicon:
            len_scaffold_to_use += 1
            scaffold_dna_search = scaffold_seq[-len_scaffold_to_use:] + extension_seq
        return scaffold_start_loc, scaffold_seq[-len_scaffold_to_use:]

    else:
        extension_seq = reverse_complement(extension_seq)
        scaffold_seq = reverse_complement(scaffold_seq)

        scaffold_start_loc = amplicon.index(extension_seq) + len(extension_seq)
        scaffold_dna_search = extension_seq + scaffold_seq[0:len_scaffold_to_use]
        while scaffold_dna_search in amplicon:
            len_scaffold_to_use += 1
            scaffold_dna_search = extension_seq + scaffold_seq[0:len_scaffold_to_use]
        return scaffold_start_loc, scaffold_seq[0:len_scaffold_to_use]


def main():
    parser = argparse.ArgumentParser(description="Summarize scaffold read through")
    parser.add_argument("-f", "--CRISPResso2_folder", type=str, help="CRISPResso output folder containing finished analysis", required=True)
    parser.add_argument("-o", "--output_folder", type=str, help="Output folder", required=True)
    parser.add_argument("-s", type=str, default="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC")
    parser.add_argument("-g", help='AA id')
    parser.add_argument("-l", type=int)

    args = parser.parse_args()

    # get AA info
    cur.execute("select sp1.bases, sp2.bases, beacon1.bases, beacon2.bases from atg_atg "
                "join modified_rna as m1 on m1.id = atg_atg.atg1 "
                "join modified_rna as m2 on m2.id = atg_atg.atg2 "
                "join atgrna as a1 on a1.id = m1.rna "
                "join atgrna as a2 on a2.id = m2.rna "
                "join dna_oligo as sp1 on sp1.id = a1.spacer "
                "join dna_oligo as sp2 on sp2.id = a2.spacer "
                "join dna_sequence as beacon1 on beacon1.id = a1.beacon "
                "join dna_sequence as beacon2 on beacon2.id = a2.beacon "
                "join registry_entity as spp on spp.id = atg_atg.spacer_pair "
                "where atg_atg.file_registry_id$ = %s", [args.g])
    sp1_seq, sp2_seq, beacon1_seq, beacon2_seq = cur.fetchone()

    #
    cs2_info = CRISPRessoShared.load_crispresso_info(args.CRISPResso2_folder)
    refs = cs2_info['results']['refs']
    stats = {}

    #  read Allele file
    if not cs2_info["running_info"]["args"].write_detailed_allele_table:
        raise Exception('CRISPResso run must be run with the parameter --write_detailed_allele_table')

    z = zipfile.ZipFile(os.path.join(args.CRISPResso2_folder, cs2_info["running_info"]["allele_frequency_table_zip_filename"]))
    zf = z.open(cs2_info["running_info"]["allele_frequency_table_filename"])
    df_alleles = pd.read_csv(zf, sep="\t")
    df_ref = df_alleles[df_alleles["Reference_Name"] == "Beacon"]
    df_ref["ref_positions"] = df_ref["ref_positions"].apply(arrstr_to_arr)

    # check scaffold read through for atg1
    scaffold_info = get_scaffold_search(refs['Beacon']['sequence'], beacon1_seq, args.s, args.l)
    print(scaffold_info)
    df_ref["scaffold_loc"] = df_ref.apply(lambda row: row["ref_positions"].index(scaffold_info[0] - 1) + 1, axis=1, result_type='expand')
    df = df_ref[df_ref.apply(lambda row: row["Aligned_Sequence"][row["scaffold_loc"]:(row["scaffold_loc"] + len(scaffold_info[1]))] == scaffold_info[1], axis=1)]
    stats["beacon1"] = df["#Reads"].sum()
    print(df)
    df.to_csv("test.txt")

    # check scaffold read through for atg2
    scaffold_info = get_scaffold_search(refs['Beacon']['sequence'], beacon2_seq, args.s, args.l)
    print(scaffold_info)
    df_ref["scaffold_loc"] = df_ref.apply(lambda row: row["ref_positions"].index(scaffold_info[0] + 1) - 1, axis=1, result_type='expand')
    print(df_ref.apply(lambda row: row["Aligned_Sequence"][(row["scaffold_loc"] - len(beacon2_seq) - len(scaffold_info[1])):row["scaffold_loc"]], axis=1, result_type='expand').iloc[19])
    print(df_ref["ref_positions"].iloc[19])
    print(df_ref["scaffold_loc"].iloc[19])
    print(reverse_complement(df_ref["Aligned_Sequence"].iloc[19]))
    df = df_ref[df_ref.apply(lambda row: row["Aligned_Sequence"][(row["scaffold_loc"] - len(beacon2_seq) - len(scaffold_info[1])):(row["scaffold_loc"] - len(beacon2_seq))] == scaffold_info[1], axis=1)]
    stats["beacon2"] = df["#Reads"].sum()

    pd.Series(stats).to_csv(os.path.join(args.CRISPResso2_folder, "CRISPResso_scaffold.stats.txt"), sep="\t", header=False, index=True, na_rep=0)


if __name__ == "__main__":
    main()
