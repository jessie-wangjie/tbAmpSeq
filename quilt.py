import argparse
import glob
import os
import pandas as pd
import quilt3
import altair as alt


def chart_heatmap(df, opts):
    file = CSV.split(".")[0]
    title = f'{opts.z}_{file}'
    return alt.Chart(df, title=title).mark_circle(size=100).encode(x=opts.y, y=opts.x,
                                                                   color=alt.Color(opts.z, scale=alt.Scale(scheme=opts.color)),
                                                                   tooltip=list(df.columns)).properties(width=550)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate quilt package', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", help='TB id')

    args = parser.parse_args()
    tbid = args.m

    # plate plot
    files = glob.glob(tbid + "/*/CRISPResso_stats.json")
    stats = pd.DataFrame()
    for s in files:
        stats = pd.concat([stats, pd.read_json(s,orient="index").T])
    print(stats)

    # Create test directories
    TEST_DIR = "test_workflow"
    SUB_DIR = "subdir"

    # create test data
#    p = quilt3.Package()
#    p.set("data.csv", "s3://tb-ngs-quilt/CRISPResso_on_9/CRISPResso_quantification_of_editing_frequency.txt",
#      meta={"type": "csv"})

#    top_hash = p.build("jwang/test_data")

#    p.push(
#        "jwang/test_data",
#        "s3://tb-ngs-quilt",
#        message="Updated version my package"
#    )

