import argparse
import glob

import altair as alt
import pandas as pd
import quilt3

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate quilt package', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", help='TB id')

    args = parser.parse_args()
    tbid = args.m

    # plate plot
    files = glob.glob(tbid + "/*/CRISPResso_stats.json")

    data = pd.DataFrame()
    for s in files:
        stats = pd.concat([data, pd.read_json(s, orient="index").T])
    data["x"] = data["well"].str.get(-1)
    data["y"] = data["well"].str.get(0)

    # plots
    selector = alt.selection_single(empty='all', fields=['sample_name'])
    base = alt.Chart(data).properties(width=250, height=250).add_selection(selector)

    # draw plate plots
    heatmap = base.mark_circle().encode(x='x:O', y='y:O', size='beacon_placement_percentage:Q',
                                        color=alt.condition(selector, 'sample_name:O', alt.value('lightgray'), legend=None))
    base.save(tbid + ".plate.json")

    # draw bar plots
    bar = base.transform_fold(["beacon_placement_percentage", "perfect_beacon_percent"]).mark_bar().encode(
        x='sample_name:O', y='value:Q', color='key:N').transform_filter(selector)

    chart = heatmap | bar
    chart.save(tbid + "report.json")

    # Create test directories
    TEST_DIR = "test_workflow"
    SUB_DIR = "subdir"

    # create test data
#    p = quilt3.Package()
#    p.set("data.csv", "s3://tb-ngs-quilt/CRISPResso_on_9/CRISPResso_quantification_of_editing_frequency.txt",
#      meta={"type": "csv"})
    #

    # edit a preexisting package
    quilt3.Package.install(
        "jwang/test_data",
        "s3://tb-ngs-quilt",
    )
    p = quilt3.Package.browse('jwang/test_data')

    # adding data
    p.set(tbid + ".report.json")

    # Saving a package manifest locally
    top_hash = p.build("jwang/test_data")
    # Pushing a package to a remote registry
    p.push(
        "jwang/test_data",
        "s3://tb-ngs-quilt",
        message="Updated version my package"
    )

