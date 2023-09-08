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


def platemap(data):
    # prepare the data for plotting
    cols = {
        "AA": ["plate", "x", "y", "samplename", "aaanid", "ppid", "total_aligned_read_num", "aligned_percentage", "beacon_placement_percentage"],
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
    d = pd.DataFrame(data)

    # draw plate plots
    bp = alt.Chart(d).mark_circle(size=300, filled=True).properties(width=300, height=200).encode(
        x=alt.X("x:Q").axis(title="").scale(domain=[1, 12]),
        y=alt.Y("y:O").axis(title="").scale(domain=["A", "B", "C", "D", "E", "F", "G", "H"]),
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
            "SG": ["plate", "x", "y", "total_read_num", "merged_r1r2_read_num", "wt_aligned_read_num", "indel_read_num"],
            "PN": ["plate", "x", "y", "total_read_num", "merged_r1r2_read_num", "wt_aligned_read_num", "PE_aligned_read_num", "PE_indel_read_num"]}

    d = pd.DataFrame(columns=["x", "y", "Variable", "Value"])
    for k, v in data.items():
        v = pd.DataFrame(v)[cols[k]]
        if k == "AA":
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


def main():
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
    entry_name, entry_url = cur.fetchone()

    # stats table
    data = {}
    writer = pd.ExcelWriter(os.path.join(input, input + ".stats.xlsx"), engine="xlsxwriter", engine_kwargs={"options": {"strings_to_numbers": True}})

    cols = {"AA": ["plate", "x", "y", "well", "animal_group", "samplename", "miseq_sample_name", "aaanid", "ppid", "spp_id", "total_read_num",
                   "merged_r1r2_read_num", "total_aligned_read_num", "aligned_percentage", "wt_aligned_read_num",
                   "beacon_aligned_read_num", "beacon_indel_read_num", "beacon_sub_read_num", "beacon_indel_percentage",
                   "beacon_sub_percentage", "wt_aligned_percentage", "beacon_placement_percentage", "perfect_beacon_percent",
                   "beacon_fidelity", "beacon_fidelity_1mm", "beacon_fidelity_2mm"],
            "SG": ["plate", "x", "y", "well", "animal_group", "samplename", "miseq_sample_name", "aaanid", "ppid", "total_read_num",
                   "merged_r1r2_read_num", "aligned_percentage", "wt_aligned_read_num", "indel_read_num", "sub_read_num",
                   "indel_percentage"],
            "PN": ["plate", "x", "y", "well", "animal_group", "samplename", "miseq_sample_name", "aaanid", "ppid", "spp_id", "total_read_num",
                   "merged_r1r2_read_num", "total_aligned_read_num", "aligned_percentage", "wt_aligned_read_num",
                   "PE_aligned_read_num", "Scaffold_aligned_read_num", "PE_indel_read_num", "PE_sub_read_num",
                   "PE_indel_percentage", "PE_sub_percentage", "wt_aligned_percentage", "PE_percentage"]}

    files = glob.glob(input + "/*/CRISPResso_quilt_stats.json")
    for f in files:
        d = json.load(open(f))
        d["x"] = int(d["well"][1:])
        d["y"] = d["well"][0]
        # d["plate"] = d["plate"] + " " + re.sub(".*(PRIP\d).*", "\\1", d["miseq_sample_name"])
        # d["plate"] = d["plate"] + " " + re.sub(".*(set\d).*", "\\1", d["miseq_sample_name"])
        # d["plate"] = d["plate"] + " " + re.sub(".*(Plate-\d).*", "\\1", d["miseq_sample_name"])
        type = d["aaanid"][0:2] if d["aaanid"][0:2] != "OT" else "SG"
        data.setdefault(type, []).append(d)

    for k, v in data.items():
        v = pd.DataFrame(v)[cols[k]].sort_values(["plate", "x", "y"])
        v.to_csv(os.path.join(input, "stats." + k + ".csv"), index=False)
        v.to_excel(writer, sheet_name=k, index=False, float_format="%.2f")

    # draw plate plot
    chart = platemap(data)

    # draw beacon fidelity plot
    if "AA" in data:
        chart = alt.hconcat(chart, fidelitymap(data["AA"])).resolve_scale(color="independent", shape="independent")
    chart.save(os.path.join(input, "platemap.json"))

    # draw alignment plots
    chart = barstats(data)
    chart.save(os.path.join(input, "alignment_stats.json"))

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
        quilt3.Package.install("AmpSeq/" + ngs_id, "s3://tb-ngs-quilt/")
        p = quilt3.Package.browse("AmpSeq/" + ngs_id)
    else:
        p = quilt3.Package()

    # adding data
    # input package
    p.set_dir("fastq/" + pipeline_run_id[:-1], pipeline_run_id[:-1])

    # output package
    preview = ["platemap.json", "alignment_stats.json"]
    for f in glob.glob(os.path.join(input, "stats.*.csv")):
        preview.append(os.path.basename(f))
        p.set(os.path.join(pipeline_run_id, os.path.basename(f)), f)

    p.set(pipeline_run_id + "/qw_stats.csv", input + "/" + "qw_stats.csv")
    p.set(pipeline_run_id + "/" + pipeline_run_id + ".stats.xlsx", input + "/" + input + ".stats.xlsx")
    p.set(pipeline_run_id + "/platemap.json", input + "/platemap.json")
    p.set(pipeline_run_id + "/alignment_stats.json", input + "/alignment_stats.json")
    p.set(pipeline_run_id + "/status.txt", input + "/" + "status.txt")
    p.set_dir(pipeline_run_id + "/cs2_alignment_html", input + "/cs2_alignment_html/")
    p.set_meta({"Benchling Entry": entry_name, "Benchling URL": entry_url})
    pd.Series(preview).to_json(input + "/quilt_summarize.json", orient="records")
    p.set(pipeline_run_id + "/quilt_summarize.json", input + "/quilt_summarize.json")

    # Pushing a package to a remote registry
    with Capturing() as output:
        p.push("AmpSeq/" + ngs_id, "s3://tb-ngs-quilt/", force=True)
    base_url = output[1].split()[-1]
    full_url = f"{base_url}/tree/{p.top_hash}"
    print(full_url)


if __name__ == "__main__":
    main()
