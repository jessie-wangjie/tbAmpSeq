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
    stats = pd.DataFrame()
    for s in files:
        stats = pd.concat([stats, pd.read_json(s, orient="index").T])
    stats["x"] = stats["well"].str.get(-1)
    stats["y"] = stats["well"].str.get(0)

    # draw plate plots
    chart = alt.Chart(stats).mark_circle().encode(
        x='x:O',
        y='y:O',
        size='beacon_placement_percentage:Q'
    )
    chart.save(tbid + "plate.json")

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
    p.set(tbid + "plate.json")

    # Saving a package manifest locally
    top_hash = p.build("jwang/test_data")
    # Pushing a package to a remote registry
    p.push(
        "jwang/test_data",
        "s3://tb-ngs-quilt",
        message="Updated version my package"
    )

