import argparse
import glob
import sys
from io import StringIO
import altair as alt
import pandas as pd
import quilt3
import re
import os
import json
from utils.base import *


class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio  # free up some memory
        sys.stdout = self._stdout


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate quilt package", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", help="TB id")
    parser.add_argument("-i", help="Ampseq result folder", default="./")

    args = parser.parse_args()
    pipeline_run_id = args.m
    input = args.i
    ngs_id = re.sub(".*(BTB\d+).*", "\\1", pipeline_run_id)
    cur.execute("select entry.name, entry.url from ngs_tracking "
                "join registration_origin on registration_origin.entity_id = ngs_tracking.id "
                "join entry on entry.id = registration_origin.origin_entry_id "
                "where file_registry_id$ = %s", [ngs_id])
    # entry_name, entry_url = cur.fetchone()

    # stats table
    data = {}
    writer = pd.ExcelWriter(os.path.join(input, input + ".stats.xlsx"), engine="xlsxwriter", engine_kwargs={"options": {"strings_to_numbers": True}})

    cols = {"WT": ["samplename", "aaanid", "ppid", "total_read_num", "merged_r1r2_read_num", "aligned_percentage", "wt_aligned_read_num"],
            "SP": ["samplename", "aaanid", "ppid", "total_read_num",
                   "merged_r1r2_read_num", "total_aligned_read_num", "aligned_percentage", "wt_aligned_read_num",
                   "beacon_aligned_read_num", "beacon_indel_read_num", "beacon_sub_read_num", "beacon_indel_percentage",
                   "beacon_sub_percentage", "wt_aligned_percentage", "beacon_placement_percentage", "perfect_beacon_percent",
                   "beacon_fidelity"],
            "AN": ["samplename", "aaanid", "ppid", "total_read_num",
                   "merged_r1r2_read_num", "total_aligned_read_num", "aligned_percentage", "wt_aligned_read_num",
                   "beacon_aligned_read_num", "beacon_indel_read_num", "beacon_sub_read_num", "beacon_indel_percentage",
                   "beacon_sub_percentage", "wt_aligned_percentage", "beacon_placement_percentage", "perfect_beacon_percent",
                   "beacon_fidelity"],
            "SG": ["plate", "x", "y", "well", "samplename", "miseq_sample_name", "aaanid", "ppid", "total_read_num",
                   "merged_r1r2_read_num", "aligned_percentage", "wt_aligned_read_num", "indel_read_num", "sub_read_num",
                   "indel_percentage"],
            "PN": ["plate", "x", "y", "well", "samplename", "miseq_sample_name", "aaanid", "ppid", "spp_id", "total_read_num",
                   "merged_r1r2_read_num", "total_aligned_read_num", "aligned_percentage", "wt_aligned_read_num",
                   "PE_aligned_read_num", "Scaffold_aligned_read_num", "PE_indel_read_num", "PE_sub_read_num",
                   "PE_indel_percentage", "PE_sub_percentage", "wt_aligned_percentage", "PE_percentage"]}

    files = glob.glob(input + "/*/CRISPResso_quilt_stats.json")
    for f in files:
        d = json.load(open(f))
        # d["x"] = int(d["well"][1:])
        # d["y"] = d["well"][0]
        # d["plate"] = d["plate"] + " " + re.sub(".*(PRIP\d).*", "\\1", d["miseq_sample_name"])
        if d["aaanid"] == "":
            type = "WT"
        elif d["aaanid"][0:2] != "OT":
            d["aaanid"][0:2]
        else:
            type = "SG"
        data.setdefault(type, []).append(d)

    for k, v in data.items():
        v = pd.DataFrame(v)[cols[k]]
        v.to_csv(os.path.join(input, "stats." + k + ".csv"), index=False)
        v.to_excel(writer, sheet_name=k, index=False, float_format="%.2f")

    # qw table
    files = glob.glob(input + "/*/CRISPResso_qw_stats.txt")
    qw_data = pd.DataFrame()
    for s in files:
        qw_data = pd.concat([qw_data, pd.read_csv(s, sep="\t")])
    if "WT" not in data.keys():
        qw_data = qw_data[
            ["samplename", "amplicon", "window_name", "window_region", "unmodified", "modified", "indels", "insertion", "deletion", "substitution",
            "whole_window_deletion"]]
    qw_data.to_csv(input + "/qw_stats.csv", index=False)
    qw_data.to_excel(writer, sheet_name="qw_stats", index=False, float_format="%.2f")
    writer.close()


    # check if the package existed
    if "AmpSeq/" + ngs_id in list(quilt3.list_packages("s3://tb-ngs-quilt/")):
        quilt3.Package.install("AmpSeq/" + ngs_id, "s3://tb-ngs-quilt/")
        p = quilt3.Package.browse("AmpSeq/" + ngs_id)
    else:
        p = quilt3.Package()

    # adding data
    # input package
    # p.set_dir("fastq/" + pipeline_run_id[:-1], pipeline_run_id[:-1])

    # output package
    preview = []
    for f in glob.glob(os.path.join(input, "stats.*.csv")):
        preview.append(os.path.basename(f))
        p.set(os.path.join(pipeline_run_id, os.path.basename(f)), f)

    p.set(pipeline_run_id + "/qw_stats.csv", input + "/" + "qw_stats.csv")
    p.set(pipeline_run_id + "/" + pipeline_run_id + ".stats.xlsx", input + "/" + input + ".stats.xlsx")
    p.set_dir(pipeline_run_id + "/cs2_alignment_html", input + "/cs2_alignment_html/")
    # p.set_meta({"Benchling Entry": entry_name, "Benchling URL": entry_url})
    pd.Series(preview).to_json(input + "/quilt_summarize.json", orient="records")
    p.set(pipeline_run_id + "/quilt_summarize.json", input + "/quilt_summarize.json")

    # Pushing a package to a remote registry
    with Capturing() as output:
        p.push("AmpSeq/" + ngs_id, "s3://tb-ngs-quilt/", force=True)
    base_url = output[1].split()[-1]
    full_url = f"{base_url}/tree/{p.top_hash}"
    print(full_url)
