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
    d = data["total_aligned_read_num"].groupby(data["miseq_sample_name"]).sum().rename("sample_aligned_read_num")
    data = pd.merge(data, d, on="miseq_sample_name")
    data = data.loc[data["sample_aligned_read_num"] > 1000]

    # draw plate plots
    base = alt.Chart(data, title="BP%").properties(width=550, height=400).encode(
        x=alt.X("x:Q").axis(title="").scale(domain=[0, 12]).stack(False),
        y=alt.Y("y:O").axis(title="").scale(domain=["A", "B", "C", "D", "E", "F", "G", "H"]),
        color=alt.Color("aaanid", scale=alt.Scale(domain=["AA1520", "AA2037", "AA2905"])))

    bp = base.mark_arc(stroke="#ffffff", innerRadius=10).encode(
        radius=alt.Radius("outer:Q", scale=alt.Scale(rangeMax=20)),
        radius2=alt.Radius2("inner:Q"),
        theta=alt.Theta("bp:Q")).transform_calculate(bp='datum.beacon_placement_percentage * PI / 50').transform_calculate(
        outer="{'AA1520': '10', 'AA2037': '15', 'AA2905': '20'}[datum.aaanid]").transform_calculate(
        inner="{'AA1520': '5', 'AA2037': '10', 'AA2905': '15'}[datum.aaanid]")

    chart = alt.layer(bp).facet(row=alt.Row("plate:O").title(""))
    return chart


def platemap_cargo(data):
    # prepare the data for plotting
    d = data["total_aligned_read_num"].groupby(data["miseq_sample_name"]).sum().rename("sample_aligned_read_num")
    data = pd.merge(data, d, on="miseq_sample_name")
    data = data.loc[data["sample_aligned_read_num"] > 1000]

    d = data.loc[data["aaanid"] == "AA1520", ["x", "y", "aaanid", "cargo_placement_percentage"]]
    d["aaanid"] = "PGI"
    data = pd.concat([data[["x", "y", "aaanid", "beacon_placement_percentage"]], d.rename(columns={"cargo_placement_percentage": "beacon_placement_percentage"})])

    # draw plate plots
    base = alt.Chart(data, title="BP%").properties(width=550, height=400).encode(
        x=alt.X("x:Q").axis(title="").scale(domain=[0, 12]).stack(False),
        y=alt.Y("y:O").axis(title="").scale(domain=["A", "B", "C", "D", "E", "F", "G", "H"]),
        color=alt.Color("aaanid", scale=alt.Scale(domain=["AA1520", "AA2037", "AA2905", "PGI"])))

    bp = base.mark_arc(stroke="#ffffff", innerRadius=10).encode(
        radius=alt.Radius("outer:Q", scale=alt.Scale(rangeMax=20)),
        radius2=alt.Radius2("inner:Q"),
        theta=alt.Theta("bp:Q")).transform_calculate(bp='datum.beacon_placement_percentage * PI / 50').transform_calculate(
        outer="{'AA1520': '10', 'AA2037': '15', 'AA2905': '20', 'PGI': 25}[datum.aaanid]").transform_calculate(
        inner="{'AA1520': '5', 'AA2037': '10', 'AA2905': '15', 'PGI': 20}[datum.aaanid]")

    chart = alt.layer(bp).facet(row=alt.Row("plate:O").title(""))
    return chart


def fidelitymap(data):
    # prepare the data for plotting
    d = data["total_aligned_read_num"].groupby(data["miseq_sample_name"]).sum().rename("sample_aligned_read_num")
    data = pd.merge(data, d, on="miseq_sample_name")
    data = data.loc[data["sample_aligned_read_num"] > 1000]

    # draw plate plots
    base = alt.Chart(data, title="BF%").properties(width=550, height=400).encode(
        x=alt.X("x:Q").axis(title="").scale(domain=[0, 12]).stack(False),
        y=alt.Y("y:O").axis(title="").scale(domain=["A", "B", "C", "D", "E", "F", "G", "H"]),
        color=alt.Color("aaanid", scale=alt.Scale(domain=["AA1520", "AA2037", "AA2905"])))

    bp = base.mark_arc(stroke="#ffffff", innerRadius=10).encode(
        radius=alt.Radius("outer:Q", scale=alt.Scale(rangeMax=20)),
        radius2=alt.Radius2("inner:Q"),
        theta=alt.Theta("bp:Q")).transform_calculate(bp='datum.beacon_fidelity * PI / 50').transform_calculate(
        outer="{'AA1520': '10', 'AA2037': '15', 'AA2905': '20'}[datum.aaanid]").transform_calculate(
        inner="{'AA1520': '5', 'AA2037': '10', 'AA2905': '15'}[datum.aaanid]")

    chart = alt.layer(bp).facet(row=alt.Row("plate:O").title(""))
    return chart


def fidelitymap_cargo(data):
    # prepare the data for plotting
    d = data["total_aligned_read_num"].groupby(data["miseq_sample_name"]).sum().rename("sample_aligned_read_num")
    data = pd.merge(data, d, on="miseq_sample_name")
    data = data.loc[data["sample_aligned_read_num"] > 1000]

    d = data.loc[data["aaanid"] == "AA1520", ["x", "y", "aaanid", "beacon_fidelity"]]
    d["aaanid"] = "PGI"
    data = pd.concat([data[["x", "y", "aaanid", "beacon_fidelity"]], d.rename(columns={"cargo_fidelity": "beacon_fidelity"})])

    # draw plate plots
    base = alt.Chart(data, title="BF%").properties(width=550, height=400).encode(
        x=alt.X("x:Q").axis(title="").scale(domain=[0, 12]).stack(False),
        y=alt.Y("y:O").axis(title="").scale(domain=["A", "B", "C", "D", "E", "F", "G", "H"]),
        color=alt.Color("aaanid", scale=alt.Scale(domain=["AA1520", "AA2037", "AA2905", "PGI"])))

    bp = base.mark_arc(stroke="#ffffff", innerRadius=10).encode(
        radius=alt.Radius("outer:Q", scale=alt.Scale(rangeMax=20)),
        radius2=alt.Radius2("inner:Q"),
        theta=alt.Theta("bp:Q")).transform_calculate(bp='datum.beacon_fidelity * PI / 50').transform_calculate(
        outer="{'AA1520': '10', 'AA2037': '15', 'AA2905': '20', 'PGI': 25}[datum.aaanid]").transform_calculate(
        inner="{'AA1520': '5', 'AA2037': '10', 'AA2905': '15', 'PGI': 20}[datum.aaanid]")

    chart = alt.layer(bp).facet(row=alt.Row("plate:O").title(""))
    return chart


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
    data = []
    writer = pd.ExcelWriter(os.path.join(input, input + ".stats.xlsx"), engine="xlsxwriter", engine_kwargs={"options": {"strings_to_numbers": True}})

    cargo_cols = ["plate", "x", "y", "well", "samplename", "miseq_sample_name", "aaanid", "ppid", "spp_id",
            "total_read_num", "merged_r1r2_read_num", "total_aligned_read_num", "aligned_percentage", "wt_aligned_read_num",
            "beacon_aligned_read_num", "beacon_indel_read_num", "beacon_sub_read_num",
            "cargo_aligned_read_num", "cargo_indel_read_num", "cargo_sub_read_num", "wt_aligned_percentage",
            "beacon_placement_percentage", "perfect_beacon_percent", "beacon_fidelity",
            "cargo_placement_percentage", "perfect_cargo_percent", "cargo_fidelity"]

    cols = ["plate", "x", "y", "well", "samplename", "miseq_sample_name", "aaanid", "ppid", "spp_id",
            "total_read_num", "merged_r1r2_read_num", "total_aligned_read_num", "aligned_percentage", "wt_aligned_read_num",
            "beacon_aligned_read_num", "beacon_indel_read_num", "beacon_sub_read_num", "wt_aligned_percentage",
            "beacon_placement_percentage", "perfect_beacon_percent", "beacon_fidelity"]

    files = glob.glob(input + "/*/CRISPResso_quilt_stats.json")
    for f in files:
        d = json.load(open(f))
        d["x"] = int(d["well"][1:])
        d["y"] = d["well"][0]
        if d["plate"] is None:
            d["plate"] = ""
        data.append(d)

    if "cargo_aligned_read_num" in data[0].keys():
        data = pd.DataFrame(data)[cargo_cols].sort_values(["plate", "x", "y", "aaanid"])
    else:
        data = pd.DataFrame(data)[cols].sort_values(["plate", "x", "y", "aaanid"])
    data.to_csv(os.path.join(input, "stats.csv"), index=False)
    data.to_excel(writer, sheet_name="AA", index=False, float_format="%.2f")

    # draw plate plot
    if "cargo_aligned_read_num" in data:
        chart = platemap_cargo(data)
        chart = alt.hconcat(chart, fidelitymap_cargo(data)).resolve_scale(color="independent", shape="independent")
    else:
        chart = platemap(data)
        chart = alt.hconcat(chart, fidelitymap(data)).resolve_scale(color="independent", shape="independent")

    # draw beacon fidelity plot

    chart.save(os.path.join(input, "platemap.json"))

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
    p.set_dir("fastq/" + pipeline_run_id[:-1], pipeline_run_id[:-1])

    # output package
    preview = ["platemap.json"]
    for f in glob.glob(os.path.join(input, "stats.csv")):
        preview.append(os.path.basename(f))
        p.set(os.path.join(pipeline_run_id, os.path.basename(f)), f)

    p.set(pipeline_run_id + "/qw_stats.csv", input + "/" + "qw_stats.csv")
    p.set(pipeline_run_id + "/" + pipeline_run_id + ".stats.xlsx", input + "/" + input + ".stats.xlsx")
    p.set(pipeline_run_id + "/platemap.json", input + "/platemap.json")
    p.set(pipeline_run_id + "/status.txt", input + "/" + "status.txt")
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

