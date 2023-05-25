import argparse
import glob
import sys
from io import StringIO
import altair as alt
import pandas as pd
import quilt3
import re
import os

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
                                                                                                     tooltip=['aaanid','ppid','spp_id',alt.Tooltip('wt_aligned_percentage',title='No beacon %'),alt.Tooltip('total_read_num',title='Total reads'),alt.Tooltip('samplename',title='Sample name')],
                                                                                                     text=alt.Text(
                                                                                                         'beacon_placement_percentage:Q',
                                                                                                         format='.0f'))

    pbp = alt.Chart(data).mark_circle(size=300).properties(width=300, height=200).encode(
        x=alt.X('x:Q', axis=alt.Axis(title=''), scale=alt.Scale(domain=[1, 12])), y=alt.Y('y:O', axis=alt.Axis(title='')),
        color=alt.Color('perfect_beacon_percent:Q', scale=alt.Scale(scheme="oranges", domain=[0, 100]),
                        legend=alt.Legend(title="Perfect BP %")))

    pbp_text = alt.Chart(data).mark_text(size=8, dx=0, dy=0, color='black', fontWeight="bold").encode(x=alt.X('x:Q'),
                                                                                                      y=alt.Y('y:O'),
                                                                                                      tooltip=['aaanid','ppid','spp_id',alt.Tooltip('wt_aligned_percentage',title='No beacon %'),alt.Tooltip('total_read_num',title='Total reads'),alt.Tooltip('samplename',title='Sample name')],
                                                                                                      text=alt.Text(
                                                                                                          'perfect_beacon_percent:Q',
                                                                                                          format='.0f'))
    chart = alt.hconcat(alt.layer(bp, bp_text).facet(row='plate:O'), alt.layer(pbp, pbp_text).facet(row='plate:O')).resolve_scale(
        color='independent')
    return chart


def barstats(data):
    # Add new calculated columns to the data
    data['complete_beacon_num'] = data['beacon_aligned_read_num'] * data['perfect_beacon_percent'] / 100
    data['imperfect_beacon_num'] = data['beacon_aligned_read_num'] - data['perfect_beacon_num']
    data['no_beacon_num'] = data['merged_r1r2_read_num'] - data['complete_beacon_num'] - data['imperfect_beacon_num']

    # Convert the data from wide to long format
    data_melted = pd.melt(data, id_vars=['x', 'y'], 
                      value_vars=['total_read_num','merged_r1r2_read_num','complete_beacon_num', 'imperfect_beacon_num', 'no_beacon_num'], 
                      var_name='Variable', value_name='Value')
    
    # Add a new column 'Stacked_Variable' for stacking 'complete_beacon_num', 'imperfect_beacon_num' and 'perfect_beacon_num' together
    data_melted['Stacked_Variable'] = data_melted['Variable'].replace({
        'no_beacon_num': 'aligned_read_num', 
        'imperfect_beacon_num': 'aligned_read_num',
        'complete_beacon_num': 'aligned_read_num'
    })

    # Define mapping from 'Variable' to legend labels
    legend_labels = {
        "total_read_num": "Total reads",
        "merged_r1r2_read_num": "Total merged R1/R2 reads",
        "complete_beacon_num": "Complete beacon insertion",
        "imperfect_beacon_num": "Imperfect beacon insertion",
        "no_beacon_num": "No beacon insertion"
    }
    data_melted['Legend'] = data_melted['Variable'].map(legend_labels)

    print(data_melted.iloc[20:40])
    # Define the base chart with common elements
    base = alt.Chart(data_melted).encode(
        alt.X('Stacked_Variable:O', title='', axis=None, 
              sort=['total_read_num', 'merged_r1r2_read_num', 'aligned_read_num']),  # Remove x-axis and sort bars
        alt.Y('Value:Q', title=''),
        alt.Color('Legend:N', legend=alt.Legend(title="Variable"))  # Add color legend
    ).properties(
        width=100,
        height=100
    )

    # Define the bar chart
    bar_chart = base.mark_bar().encode(
        color='Legend:N'
    )

    # Arrange the charts in a grid based on 'x' and 'y'
    final_chart = bar_chart.facet(
        column='x:N',
        row='y:N'
    )

    return final_chart


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate quilt package', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", help='TB id')
    parser.add_argument("-i", help='Ampseq result folder', default="./")

    args = parser.parse_args()
    pipeline_run_id = args.m
    input = args.i
    ngs_id = re.sub(".*(BTB\d+).*", "\\1", pipeline_run_id)

    # stats table
    writer = pd.ExcelWriter(input + "/stats.xlsx", engine="xlsxwriter", engine_kwargs={'options': {'strings_to_numbers': True}})

    files = glob.glob(input + "/*/CRISPResso_quilt_stats.json")
    data = pd.DataFrame()
    for s in files:
        data = pd.concat([data, pd.read_json(s, orient="index").T])
    data["x"] = data["well"].str.extract(r"(\d+)")
    data["x"] = data["x"].astype('int')
    data["y"] = data["well"].str.get(0)
    data["plate"] = data["plate"].fillna("Plate 1")
    if "beacon_placement_percentage" in data:
        data = data[
            ["plate", "x", "y", "well", "samplename", "miseq_sample_name", "aaanid", "ppid", "spp_id", "total_read_num", "merged_r1r2_read_num",
             "total_aligned_read_num", "aligned_percentage", "wt_aligned_read_num", "beacon_aligned_read_num", "beacon_indel_read_num",
             "beacon_sub_read_num", "beacon_indel_percentage", "beacon_sub_percentage", "wt_aligned_percentage", "beacon_placement_percentage",
             "perfect_beacon_percent", "beacon_fidelity"]]
        data.to_csv(input + "/stats.csv", index=False)
        data.to_excel(writer, sheet_name="AA_AN", index=False, float_format="%.2f")
    elif "indel_percentage" in data:
        data = data[
            ["plate", "x", "y", "well", "samplename", "miseq_sample_name", "aaanid", "ppid", "total_read_num", "merged_r1r2_read_num",
             "aligned_percentage", "wt_aligned_read_num", "indel_read_num", "sub_read_num", "indel_percentage"]]
        data.to_csv(input + "/stats.sg.csv", index=False)
        data.to_excel(writer, sheet_name="SG", index=False, float_format="%.2f")
    elif "PE_percentage" in data:
        data = data[
            ["plate", "x", "y", "well", "samplename", "miseq_sample_name", "aaanid", "ppid", "spp_id", "total_read_num", "merged_r1r2_read_num",
             "total_aligned_read_num", "aligned_percentage", "wt_aligned_read_num", "PE_aligned_read_num", "Scaffold_aligned_read_num",
             "PE_indel_read_num", "PE_sub_read_num", "PE_indel_percentage", "PE_sub_percentage", "wt_aligned_percentage", "PE_percentage"]]
        data.to_csv(input + "/stats.pe.csv", index=False)
        data.to_excel(writer, sheet_name="PE", index=False, float_format="%.2f")


    # plots
    # draw plate plots
    chart = platemap(data)
    chart.save(input + "/platemap.json")

    # draw alignment plots
    chart = barstats(data)
    chart.save(input + "/alignment_stats.json")

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
    if os.path.exists(os.path.join(pipeline_run_id,"stats.csv")):
        p.set(pipeline_run_id + "/stats.csv", input + "/stats.csv")
    if os.path.exists(os.path.join(pipeline_run_id,"stats.sg.csv")):
        p.set(pipeline_run_id + "/stats.sg.csv", input + "/stats.sg.csv")
    if os.path.exists(os.path.join(pipeline_run_id,"stats.pe.csv")):
        p.set(pipeline_run_id + "/stats.pe.csv", input + "/stats.pe.csv")
    p.set(pipeline_run_id + "/qw_stats.csv", input + "/" + "qw_stats.csv")
    p.set(pipeline_run_id + "/stats.xlsx", input + "/" + "stats.xlsx")
    p.set(pipeline_run_id + "/platemap.json", input + "/platemap.json")
    p.set(pipeline_run_id + "/alignment_stats.json", input + "/alignment_stats.json")
    p.set(pipeline_run_id + "/status.txt", input + "/" + "status.txt")
    p.set_dir(pipeline_run_id + "/cs2_alignment_html", input + "/cs2_alignment_html/")
    preview = pd.Series(["status.txt", "platemap.json", "stats.csv"])
    preview.to_json(input + "/quilt_summarize.json", orient="records")
    p.set(pipeline_run_id + "/quilt_summarize.json", input + "/quilt_summarize.json")

    # Pushing a package to a remote registry
    with Capturing() as output:
        p.push("AmpSeq/" + ngs_id, "s3://tb-ngs-quilt/", force=True)
    base_url = output[1].split()[-1]
    full_url = f"{base_url}/tree/{p.top_hash}"
    print(full_url)
