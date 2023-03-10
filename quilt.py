import argparse
import glob
import sys
from io import StringIO

import altair as alt
import pandas as pd
import quilt3


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

    bp_text = alt.Chart(data).mark_text(size=8, dx=0, dy=0, color='black', fontWeight="bold").encode(
        x=alt.X('x:Q'),
        y=alt.Y('y:O'),
        text=alt.Text('beacon_placement_percentage:Q', format='.0f'))

    pbp = alt.Chart(data).mark_circle(size=300).properties(width=300, height=200).encode(
        x=alt.X('x:Q', axis=alt.Axis(title=''), scale=alt.Scale(domain=[1, 12])),
        y=alt.Y('y:O', axis=alt.Axis(title='')),
        color=alt.Color('perfect_beacon_percent:Q', scale=alt.Scale(scheme="oranges", domain=[0, 100]),
                        legend=alt.Legend(title="Perfect BP %")))

    pbp_text = alt.Chart(data).mark_text(size=8, dx=0, dy=0, color='black', fontWeight="bold").encode(
        x=alt.X('x:Q'),
        y=alt.Y('y:O'),
        text=alt.Text('perfect_beacon_percent:Q', format='.0f'))

    chart = alt.hconcat(alt.layer(bp, bp_text).facet(row='plate:O'),
                        alt.layer(pbp, pbp_text).facet(row='plate:O')).resolve_scale(color='independent')
    return chart


def barstats(data):
    # draw alignment stats
    bar = alt.Chart(data, width=60, height=alt.Step(8)).transform_fold(
        ["total_read_num", "merged_r1r2_read_num"], as_=["key", "num"]).mark_bar().encode(
        x=alt.X('num:Q', axis=alt.Axis(title='')),
        y=alt.Y('key:N', axis=None, scale=alt.Scale(domain=["total_read_num", "merged_r1r2_read_num"])),
        color=alt.Color('key:N', sort=["total_read_num", "merged_r1r2_read_num"], legend=alt.Legend(title="")))

    bar2 = alt.Chart(data).transform_fold(
        ["wt_aligned_read_num", "beacon_aligned_read_num"], as_=["key", "num"]).transform_stack(
        stack="num", as_=["x1", "x2"], groupby=["x", "y"]).mark_bar().encode(
        x=alt.X('x1:Q'),
        x2="x2:Q",
        y=alt.Y("plate:O", scale=alt.Scale(domain=["Plate 1"])),
        color=alt.Color('key:N', sort=["total_read_num", "merged_r1r2_read_num", "beacon_aligned_read_num",
                                       "wt_aligned_read_num"]))
    chart = alt.layer(bar + bar2).facet(row=alt.Row('y:O', title=""), column=alt.Column('x:O', title=""))
    return chart


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate quilt package',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", help='TB id')
    parser.add_argument("-i", help='Ampseq result folder', default="./")

    args = parser.parse_args()
    tbid = args.m
    input = args.i

    # plate plot
    files = glob.glob(input + "/*/CRISPResso_stats.json")

    data = pd.DataFrame()
    for s in files:
        data = pd.concat([data, pd.read_json(s, orient="index").T])
    data["x"] = data["well"].str.extract(r"(\d+)")
    data["x"] = data["x"].astype('int')
    data["y"] = data["well"].str.get(0)
    data["plate"] = data["plate"].fillna("Plate 1")
    data = data[["plate", "x", "y", "well", "samplename", "miseq_sample_name", "aaanid", "ppid", "total_read_num",
                 "merged_r1r2_read_num", "aligned_percentage", "wt_aligned_read_num", "beacon_aligned_read_num",
                 "wt_aligned_percentage", "beacon_placement_percentage", "beacon_indel_read_num", "beacon_sub_read_num",
                 "beacon_indel_percentage", "beacon_sub_percentage", "perfect_beacon_percent"]]
    data.to_csv(input + "/stats.csv", index=False)

    # plots
    # draw plate plots
    chart = platemap(data)
    chart.save(input + "/platemap.json")

    # draw alignment plots
    chart = barstats(data)
    chart.save(input + "/alignment_stats.json")

    # push data to quilt
    p = quilt3.Package()

    # edit a preexisting package
    #    quilt3.Package.install(
    #        "jwang/" + tbid,
    #        "s3://tb-ngs-quilt",
    #    )
    #    p = quilt3.Package.browse("jwang/" + tbid)

    # adding data
    p.set("stats.csv", input + "/stats.csv")
    p.set("platemap.json", input + "/platemap.json")
    p.set("alignment_stats.json", input + "/alignment_stats.json")
    p.set("status.txt", input + "/" + tbid + ".status.txt")
    p.set_dir("cs2_alignment_html", input + "/cs2_alignment_html/")
    preview = pd.Series(["platemap.json", "alignment_stats.json", "status.txt", "stats.csv"])
    preview.to_json(input + "/quilt_summarize.json", orient="records")
    p.set("quilt_summarize.json", input + "/quilt_summarize.json")

    # Pushing a package to a remote registry
    with Capturing() as output:
        p.push(
            "jwang/" + tbid,
            "s3://tb-quilt-test",
        )
    base_url = output[1].split()[-1]
    full_url = f"{base_url}/tree/{p.top_hash}"
    print(full_url)
