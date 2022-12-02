import argparse
import glob

import altair as alt
import pandas as pd
import quilt3

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
    data["x"] = data["well"].str.get(-1)
    data["y"] = data["well"].str.get(0)
    data.to_csv(tbid + "/stats.csv")

    # plots
    selector = alt.selection_single(empty='all', fields=['sample_name'])
    base = alt.Chart(data)

    # draw plate plots
    heatmap = base.mark_circle().properties(width=250, height=250).encode(x='x:O', y='y:O',
                size='beacon_placement_percentage:Q',
                color=alt.condition(selector, 'sample_name:O',alt.value('lightgray'), legend=None)).add_selection(selector)

    # draw bar plots
    bar = base.transform_fold(["beacon_placement_percentage", "perfect_beacon_percent"],
                              as_=['beacon', 'percent']).mark_bar(opacity=0.7).encode(
        x='sample_name:O', y=alt.Y('percent:Q', stack=None),
        color=alt.condition(selector, alt.Color("beacon:N"), alt.ColorValue("grey"))).add_selection(selector)

    chart = alt.hconcat(heatmap, bar).resolve_scale(color="independent")
    chart.save(tbid + "/report.json")

    # create test data
    p = quilt3.Package()
#    p.set("data.csv", "s3://tb-ngs-quilt/CRISPResso_on_9/CRISPResso_quantification_of_editing_frequency.txt", meta={"type": "csv"})

    # edit a preexisting package
#    quilt3.Package.install(
#        "jwang/" + tbid,
#        "s3://tb-ngs-quilt",
#    )
#    p = quilt3.Package.browse("jwang/" + tbid)

    # adding data
    p.set(tbid + "/stats.csv")
    p.set(tbid + "/report.json")
    p.set_dir(tbid + "/alignment_html")

    # Pushing a package to a remote registry
#    p.push(
#        "jwang/" + tbid,
#        "s3://tb-ngs-quilt",
#        message="Updated version my package"
#    )
