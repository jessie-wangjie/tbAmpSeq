import argparse
import glob
import sys
from io import StringIO

import altair as alt
import pandas as pd


class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate quilt package',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", help='TB id')

    args = parser.parse_args()
    tbid = args.m

    # plate plot
    files = glob.glob(tbid + "/*/CRISPResso_stats.json")

    data = pd.DataFrame()
    for s in files:
        data = pd.concat([data, pd.read_json(s, orient="index").T])
    data["x"] = data["well"].str.extract(r"(\d+)")
    data["x"] = data["x"].astype('int')
    data["y"] = data["well"].str.get(0)
    data = data[["plate", "x", "y", "well", "sample_name", "miseq_sample_name", "aaan_id", "ppid", "total_read_num",
                 "merged_r1r2_read_num", "aligned_percentage", "wt_aligned_read_num", "beacon_aligned_read_num",
                 "wt_aligned_percentage", "beacon_placement_percentage", "beacon_indel_read_num", "beacon_sub_read_num",
                 "beacon_indel_percentage", "beacon_sub_percentage", "perfect_beacon_percent"]]
    data.to_csv(tbid + "/stats.csv", index=False)

    # plots
    selector = alt.selection_single(empty='all', fields=['sample_name'])
    base = alt.Chart(data)

    # draw plate plots
    heatmap = base.mark_circle().properties(width=300, height=200).encode(x=alt.X('x:O'), y='y:O',
                size='beacon_placement_percentage:Q')

    heatmap2 = heatmap.mark_point().encode(size='perfect_beacon_percent:Q',
                                           color=alt.ColorValue("orange"))

    chart = alt.hconcat(heatmap+heatmap2,)
    chart.save(tbid + "/report.json")

    # draw bar plots
#    bar = base.transform_fold(["beacon_placement_percentage", "perfect_beacon_percent"],
#                              as_=['beacon', 'percent']).mark_bar(opacity=0.7).properties(width=600, height=200).encode(
#        x='sample_name:O', y=alt.Y('percent:Q', stack=None),
#        color=alt.condition(selector, alt.Color("beacon:N"), alt.ColorValue("grey"))).add_selection(selector)

#    chart = alt.hconcat(heatmap, bar).resolve_scale(color="independent")
#    chart.save(tbid + "/report.json")

    # create test data
#    p = quilt3.Package()
#    p.set("data.csv", "s3://tb-ngs-quilt/CRISPResso_on_9/CRISPResso_quantification_of_editing_frequency.txt", meta={"type": "csv"})

    # edit a preexisting package
#    quilt3.Package.install(
#        "jwang/" + tbid,
#        "s3://tb-ngs-quilt",
#    )
#    p = quilt3.Package.browse("jwang/" + tbid)

    # adding data
#    p.set("stats.csv", tbid + "/stats.csv")
#    p.set("report.json", tbid + "/report.json")
#    p.set_dir("alignment_html", tbid + "/alignment_html/")
#    preview = pd.Series(["report.json", "stats.csv"])
#    preview.to_json(tbid + "/quilt_summarize.json", orient="records")
#    p.set("quilt_summarize.json", tbid + "/quilt_summarize.json")

    # Pushing a package to a remote registry
#    with Capturing() as output:
#        p.push(
#            "jwang/" + tbid,
#            "s3://tb-ngs-quilt",
#        )
#    base_url = output[1].split()[-1]
#    full_url = f"{base_url}/tree/{p.top_hash}"
#    print(full_url)