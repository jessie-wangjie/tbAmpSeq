import argparse
import glob
import sys
from io import StringIO
import altair as alt
import pandas as pd
import quilt3
import re

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
    # draw plate plots
    bp = alt.Chart(data).mark_circle(size=300).properties(width=300, height=200).encode(
        x=alt.X('x:Q', axis=alt.Axis(title=''), scale=alt.Scale(domain=[1, 12])),
        y=alt.Y('y:O', axis=alt.Axis(title='')),
        color=alt.Color('beacon_placement_percentage:Q', scale=alt.Scale(scheme="blues", domain=[0, 100]),
                        legend=alt.Legend(title="BP %")))

    bp_text = alt.Chart(data).mark_text(size=8, dx=0, dy=0, color='black', fontWeight="bold").encode(x=alt.X('x:Q'),
                                                                                                     y=alt.Y('y:O'),
                                                                                                     text=alt.Text(
                                                                                                         'beacon_placement_percentage:Q',
                                                                                                         format='.0f'))

    pbp = alt.Chart(data).mark_circle(size=300).properties(width=300, height=200).encode(
        x=alt.X('x:Q', axis=alt.Axis(title=''), scale=alt.Scale(domain=[1, 12])), y=alt.Y('y:O', axis=alt.Axis(title='')),
        color=alt.Color('perfect_beacon_percent:Q', scale=alt.Scale(scheme="oranges", domain=[0, 100]),
                        legend=alt.Legend(title="Perfect BP %")))

    pbp_text = alt.Chart(data).mark_text(size=8, dx=0, dy=0, color='black', fontWeight="bold").encode(x=alt.X('x:Q'),
                                                                                                      y=alt.Y('y:O'),
                                                                                                      text=alt.Text(
                                                                                                          'perfect_beacon_percent:Q',
                                                                                                          format='.0f'))

    chart = alt.hconcat(alt.layer(bp, bp_text).facet(row='plate:O'), alt.layer(pbp, pbp_text).facet(row='plate:O')).resolve_scale(
        color='independent')
    return chart


def barstats(data):
    # draw alignment stats
    bar = alt.Chart(data, width=60, height=alt.Step(8)).transform_fold(["total_read_num", "merged_r1r2_read_num"],
                                                                       as_=["key", "num"]).mark_bar().encode(
        x=alt.X('num:Q', axis=alt.Axis(title='')),
        y=alt.Y('key:N', axis=None, scale=alt.Scale(domain=["total_read_num", "merged_r1r2_read_num"])),
        color=alt.Color('key:N', sort=["total_read_num", "merged_r1r2_read_num"], legend=alt.Legend(title="")))

    bar2 = alt.Chart(data).transform_fold(["wt_aligned_read_num", "beacon_aligned_read_num"], as_=["key", "num"]).transform_stack(
        stack="num", as_=["x1", "x2"], groupby=["x", "y"]).mark_bar().encode(x=alt.X('x1:Q'), x2="x2:Q", y=alt.Y("plate:O",
                                                                                                                 scale=alt.Scale(
                                                                                                                     domain=[
                                                                                                                         "Plate 1"])),
                                                                             color=alt.Color('key:N', sort=["total_read_num",
                                                                                                            "merged_r1r2_read_num",
                                                                                                            "beacon_aligned_read_num",
                                                                                                            "wt_aligned_read_num"]))
    chart = alt.layer(bar + bar2).facet(row=alt.Row('y:O', title=""), column=alt.Column('x:O', title=""))
    return chart


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate quilt package', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", help='TB id')
    parser.add_argument("-i", help='Ampseq result folder', default="./")

    args = parser.parse_args()
    pipeline_run_id = args.m
    input = args.i
    ngs_id = re.sub(".*(BTB\d+).*", "\\1", pipeline_run_id)

    # plate plot
    files = glob.glob(input + "/*/CRISPResso_quilt_stats.json")

    data = pd.DataFrame()
    for s in files:
        data = pd.concat([data, pd.read_json(s, orient="index").T])
    print(data)
    data["x"] = data["well"].str.extract(r"(\d+)")
    data["x"] = data["x"].astype('int')
    data["y"] = data["well"].str.get(0)
    data["plate"] = data["plate"].fillna("Plate 1")
    data = data[["plate", "x", "y", "well", "samplename", "miseq_sample_name", "aaanid", "ppid", "spp_id", "total_read_num", "merged_r1r2_read_num",
                 "aligned_percentage", "wt_aligned_read_num", "beacon_aligned_read_num", "beacon_indel_read_num", "beacon_sub_read_num",
                 "beacon_indel_percentage", "beacon_sub_percentage", "wt_aligned_percentage", "beacon_placement_percentage",
                 "perfect_beacon_percent", "beacon_fidelity"]]
    data.to_csv(input + "/stats.csv", index=False)

    # plots
    # draw plate plots
    chart = platemap(data)
    chart.save(input + "/platemap.json")

    # draw alignment plots
    chart = barstats(data)
    chart.save(input + "/alignment_stats.json")

    # check if the package existed
    if "jwang/" + ngs_id in list(quilt3.list_packages("s3://tb-quilt-test/")):
        quilt3.Package.install("jwang/" + ngs_id, "s3://tb-quilt-test/")
        p = quilt3.Package.browse("jwang/" + ngs_id)
    else:
        p = quilt3.Package()

    # adding data
    # input package
    p.set_dir("fastq/" + pipeline_run_id[:-1], pipeline_run_id[:-1])

    # output package
    p.set(pipeline_run_id + "/stats.csv", input + "/stats.csv")
    p.set(pipeline_run_id + "/platemap.json", input + "/platemap.json")
    p.set(pipeline_run_id + "/alignment_stats.json", input + "/alignment_stats.json")
    p.set(pipeline_run_id + "/status.txt", input + "/" + "status.txt")
    p.set_dir(pipeline_run_id + "/cs2_alignment_html", input + "/cs2_alignment_html/")
    preview = pd.Series(["status.txt", "platemap.json", "alignment_stats.json", "stats.csv"])
    preview.to_json(input + "/quilt_summarize.json", orient="records")
    p.set(pipeline_run_id + "/quilt_summarize.json", input + "/quilt_summarize.json")

    # Pushing a package to a remote registry
    with Capturing() as output:
        p.push("jwang/" + ngs_id, "s3://tb-quilt-test/", force=True)
    base_url = output[1].split()[-1]
    full_url = f"{base_url}/tree/{p.top_hash}"
    print(full_url)
