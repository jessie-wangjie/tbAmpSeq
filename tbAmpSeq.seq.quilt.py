import argparse
import glob
import sys
from io import StringIO
import altair as alt
import pandas as pd
import quilt3
import re
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


def platemap(data):
    # prepare the data for plotting
    cols = {
        "AA": ["plate", "x", "y", "samplename", "aaanid", "ppid", "total_aligned_read_num", "aligned_percentage", "beacon_placement_percentage"],
        "LM": ["plate", "x", "y", "samplename", "aaanid", "ppid", "total_aligned_read_num", "aligned_percentage", "beacon_placement_percentage"],
        "SG": ["plate", "x", "y", "samplename", "aaanid", "ppid", "wt_aligned_read_num", "aligned_percentage", "indel_percentage"],
        "PN": ["plate", "x", "y", "samplename", "aaanid", "ppid", "total_aligned_read_num", "aligned_percentage", "PE_percentage"]}

    d = pd.DataFrame(
        columns=["type", "plate", "x", "y", "samplename", "aaanid", "ppid", "total_aligned_read_num", "aligned_percentage", "edit_percentage"])

    for k, v in data.items():
        v = pd.DataFrame(v)[cols[k]]
        v.insert(0, "type", k)
        v.columns = d.columns
        d = pd.concat([d, v], axis=0, ignore_index=True)

    # draw plate plots
    bp = alt.Chart(d).mark_point(size=300, filled=True).properties(width=300, height=200).encode(
        x=alt.X("x:Q").axis(title="").scale(domain=[1, 12]),
        y=alt.Y("y:O").axis(title="").scale(domain=["A", "B", "C", "D", "E", "F", "G", "H"]),
        shape=alt.Shape("type:N").title("Type").legend(symbolFillColor="steelblue", symbolStrokeWidth=0),
        color=alt.Color("edit_percentage:Q").scale(scheme="blues", domain=[0, 100]).title("Editing %").legend(gradientLength=95))

    text = alt.Chart(d).mark_text(size=8, dx=0, dy=0, color="black", fontWeight="bold").encode(
        x=alt.X("x:Q"),
        y=alt.Y("y:O"),
        text=alt.Text("edit_percentage:Q", format=".0f"),
        tooltip=[alt.Tooltip("samplename", title="Sample name"),
                 alt.Tooltip("aaanid", title="AA/AN/SG/PN id"),
                 alt.Tooltip("ppid", title="PP id"),
                 alt.Tooltip("total_aligned_read_num", title="Total aligned reads"),
                 alt.Tooltip("aligned_percentage", title="Aligned %")])

    chart = alt.layer(bp, text).facet(row=alt.Row("plate:O").title(""))
    return chart


def fidelitymap(data):
    # prepare the data for plotting
    d = pd.DataFrame()
    for k, v in data.items():
        if k != "AA" and k != "LM":
            continue
        v = pd.DataFrame(v)
        v.insert(0, "type", k)
        d = d.append(v)

    d.to_csv("test.csv")
    # draw plate plots
    bp = alt.Chart(d).mark_point(size=300, filled=True).properties(width=300, height=200).encode(
        x=alt.X("x:Q").axis(title="").scale(domain=[1, 12]),
        y=alt.Y("y:O").axis(title="").scale(domain=["A", "B", "C", "D", "E", "F", "G", "H"]),
        shape=alt.Shape("type:N").title("Type").legend(symbolFillColor="#ff7f0e", symbolStrokeWidth=0),
        color=alt.Color("beacon_fidelity:Q").scale(scheme="oranges", domain=[0, 100]).title("Beacon fidelity %").legend(gradientLength=95))

    text = alt.Chart(d).mark_text(size=8, dx=0, dy=0, color="black", fontWeight="bold").encode(
        x=alt.X("x:Q"),
        y=alt.Y("y:O"),
        text=alt.Text("beacon_fidelity:Q", format=".0f"),
        tooltip=[alt.Tooltip("samplename", title="Sample name"),
                 alt.Tooltip("aaanid", title="AA/AN/SG/PN id"),
                 alt.Tooltip("perfect_beacon_percent", title="Perfect BP %")])

    chart = alt.layer(bp, text).facet(row=alt.Row("plate:O").title(""))
    return chart


def barstats(data):
    cols = {"AA": ["plate", "x", "y", "total_read_num", "merged_r1r2_read_num", "wt_aligned_read_num", "beacon_aligned_read_num",
                   "beacon_indel_read_num"],
            "LM": ["plate", "x", "y", "total_read_num", "merged_r1r2_read_num", "wt_aligned_read_num", "beacon_aligned_read_num",
                   "beacon_indel_read_num"],
            "SG": ["plate", "x", "y", "total_read_num", "merged_r1r2_read_num", "wt_aligned_read_num", "indel_read_num"],
            "PN": ["plate", "x", "y", "total_read_num", "merged_r1r2_read_num", "wt_aligned_read_num", "PE_aligned_read_num", "PE_indel_read_num"]}

    d = pd.DataFrame(columns=["x", "y", "Variable", "Value"])
    for k, v in data.items():
        v = pd.DataFrame(v)[cols[k]]
        if k == "AA" or k == "LM":
            v["perfect_beacon_read_num"] = v["beacon_aligned_read_num"] - v["beacon_indel_read_num"]
            v = pd.melt(v, id_vars=["plate", "x", "y"],
                        value_vars=["total_read_num", "merged_r1r2_read_num", "wt_aligned_read_num", "perfect_beacon_read_num",
                                    "beacon_indel_read_num"], var_name="Variable", value_name="Value")
        elif k == "SG":
            v["perfect_wt_read_num"] = v["wt_aligned_read_num"] - v["indel_read_num"]
            v = pd.melt(v, id_vars=["plate", "x", "y"],
                        value_vars=["total_read_num", "merged_r1r2_read_num", "perfect_wt_read_num", "indel_read_num"], var_name="Variable",
                        value_name="Value")
        d = pd.concat([d, v], axis=0, ignore_index=True)

    d["Stacked_Variable"] = d["Variable"].replace({
        "wt_aligned_read_num": "aligned_read_num", "perfect_beacon_read_num": "aligned_read_num", "beacon_indel_read_num": "aligned_read_num",
        "perfect_wt_read_num": "aligned_read_num", "indel_read_num": "aligned_read_num"})

    legend_labels = {"total_read_num": "Total reads", "merged_r1r2_read_num": "Total merged reads",
                     "wt_aligned_read_num": "AA/AN: WT", "perfect_beacon_read_num": "AA/AN: Perfect Beacon",
                     "beacon_indel_read_num": "AA/AN: Imperfect Beacon",
                     "perfect_wt_read_num": "SG: Perfect WT", "indel_read_num": "SG: Indels"}
    d["Legend"] = d["Variable"].map(legend_labels)

    # Define the base chart with common elements
    bar = alt.Chart(d).mark_bar().properties(width=400, height=45).encode(
        alt.X("x:O").title("").scale(domain=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]),
        alt.Y("Value:Q").title(""),
        alt.XOffset("Stacked_Variable:O").sort(["total_read_num", "merged_r1r2_read_num", "aligned_read_num"]),
        alt.Color("Legend:N").scale(
            domain=["Total reads", "Total merged reads", "SG: Perfect WT", "SG: Indels", "AA/AN: WT", "AA/AN: Imperfect Beacon",
                    "AA/AN: Perfect Beacon"]))
    chart = bar.facet(row=alt.Row("y:O").title(""), column=alt.Column("plate:O").title(""), spacing=10)
    return chart


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate quilt package", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", help="TB id")
    parser.add_argument("-i", help="Ampseq result folder", default="./")

    args = parser.parse_args()
    pipeline_run_id = args.m
    input = args.i

    ngs_id = re.sub(".*(CTB\d+).*", "\\1", pipeline_run_id)

    # stats table
    data = []
    writer = pd.ExcelWriter(os.path.join(input, input + ".stats.xlsx"), engine="xlsxwriter", engine_kwargs={"options": {"strings_to_numbers": True}})

    cols = ["samplename", "aaanid", "ppid", "total_read_num", "merged_r1r2_read_num", "total_aligned_read_num", "aligned_percentage",
            "wt_aligned_read_num", "beacon_aligned_read_num", "beacon_indel_read_num", "beacon_sub_read_num", "beacon_indel_percentage",
            "beacon_sub_percentage", "wt_aligned_percentage", "beacon_placement_percentage", "perfect_beacon_percent",
            "beacon_fidelity", "beacon_fidelity_1mm", "beacon_fidelity_2mm"]

    files = glob.glob(input + "/*/CRISPResso_quilt_stats.json")

    for f in files:
        d = json.load(open(f))
        data.append(d)

    v = pd.DataFrame(data)
    v.to_csv(os.path.join(input, "stats.csv"), index=False)
    v.to_excel(writer, index=False, float_format="%.2f")

    # qw table
    files = glob.glob(input + "/*/CRISPResso_qw_stats.txt")
    qw_data = pd.DataFrame()
    for s in files:
        qw_data = pd.concat([qw_data, pd.read_csv(s, sep="\t")])
    qw_data = qw_data[
        ["samplename", "amplicon", "window_name", "window_region", "unmodified", "modified", "indels", "insertion", "deletion", "substitution",
         "whole_window_deletion"]]
    qw_data.to_csv(input + "/qw_stats.csv", index=False)
    qw_data.to_excel(writer, sheet_name="qw_stats", index=False, float_format="%.2f")
    writer.close()

    # check if the package existed
    if "AmpSeq/" + ngs_id in list(quilt3.list_packages("s3://tb-ngs-quilt/")):
        with Capturing() as tmp:
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
    pd.Series(preview).to_json(input + "/quilt_summarize.json", orient="records")
    p.set(pipeline_run_id + "/quilt_summarize.json", input + "/quilt_summarize.json")

    # Pushing a package to a remote registry
    with Capturing() as output:
        p.push("AmpSeq/" + ngs_id, "s3://tb-ngs-quilt/", force=True)
    base_url = output[1].split()[-1]
    full_url = f"{base_url}/tree/{p.top_hash}"
    print(full_url)

