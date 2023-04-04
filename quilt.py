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


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate quilt package', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    data = data[["plate", "x", "y", "well", "samplename", "miseq_sample_name", "aaanid", "ppid", "total_read_num", "merged_r1r2_read_num",
                 "aligned_percentage", "wt_aligned_read_num", "beacon_aligned_read_num", "wt_aligned_percentage", "beacon_placement_percentage",
                 "beacon_indel_read_num", "beacon_sub_read_num", "beacon_indel_percentage", "beacon_sub_percentage", "perfect_beacon_percent"]]
    data.to_csv(tbid + "/stats.csv", index=False)

    # plots
    # draw plate plots
    heatmap = alt.Chart(data).mark_rect().properties(width=300, height=200).encode(x=alt.X('x:O'), y='y:O',
                                                                                   color=alt.Color('beacon_placement_percentage:Q',
                                                                                                   scale=alt.Scale(scheme='viridis', domain=[0, 1])))

    heatmap2 = heatmap.mark_circle().encode(size=alt.Size('perfect_beacon_percent:Q', scale=alt.Scale(domain=[0, 1])),
                                            color=alt.ColorValue("#ff7f01"))

    chart = alt.hconcat(heatmap + heatmap2).resolve_scale(color='independent')
    chart.save(tbid + "/report.json")

    # create test data
    p = quilt3.Package()

    # adding data
    p.set("stats.csv", tbid + "/stats.csv")
    p.set("report.json", tbid + "/report.json")
    p.set_dir("alignment_html", tbid + "/alignment_html/")
    preview = pd.Series(["report.json", "stats.csv"])
    preview.to_json(tbid + "/quilt_summarize.json", orient="records")
    p.set("quilt_summarize.json", tbid + "/quilt_summarize.json")

    # Pushing a package to a remote registry
    with Capturing() as output:
        p.push(
            "jwang/" + tbid,
            "s3://tb-ngs-quilt",
        )
    base_url = output[1].split()[-1]
    full_url = f"{base_url}/tree/{p.top_hash}"
    print(full_url)
