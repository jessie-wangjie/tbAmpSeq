"""
Smith-Waterman Aligner Wrapper

This module provide functions to traverse input FastA/Q file and align each
read to a template sequence with Smith-Waterman aligner, which can be found
in the clib module.
"""
import argparse
import os
import pandas as pd
import zipfile
import sys
from collections import Counter
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPRessoPlot

ALIGN_TEMPLATE = "<tr><td>" \
                 "<span class=\"template\">{prefix}</span>" \
                 "<span class=\"highlight\">{seq_to_highlight}</span>" \
                 "<span class=\"template\">{postfix}</span>" \
                 "</td><td></td><td></td></tr>"
ALIGN_PADDING = "<span class=\"padding\">{seq}</span>"
ALIGN_MATCH = "{seq}"
ALIGN_MISMATCH = "<span class=\"mm{nt}\">{nt}</span>"
ALIGN_INSERTION = "<span class=\"insertion\" len=\"{len}\" seq=\"{seq}\">{pre}</span>"
ALIGN_DELETION = "<span class=\"deletion\">{seq}</span>"
ALIGN_SOFTCLIP = "<span class=\"softclip\" len=\"{len}\">{seq}</span>"
ALIGN_RECORD = "<tr><td>{alignment}</td><td>{count}</td><td>{frac:.2f}%</td></tr>"

HTML_HEADER = '''
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Alignment Output</title>
<link 
    rel="stylesheet" 
    href="https://cdn.jsdelivr.net/npm/bootstrap@4.6.0/dist/css/bootstrap.min.css" 
    integrity="sha384-B0vP5xmATw1+K9KRQjQERJvTumQW0nPEzvF6L/Z6nronJ3oUOFUFpCjEUQouq2+l" 
    crossorigin="anonymous">
<style>
    .alignment { font-family: 'PT Mono', 'Courier New', monospace; white-space: nowrap; }
    .highlight { color: red; font-weight: bold }
    .template { color: blue; background-color: #DDD; font-weight: bold; }
    .mmA { color: white; background-color: #F6511D; }
    .mmC { color: white; background-color: #FFB400; }
    .mmG { color: white; background-color: #00A6ED; }
    .mmT { color: white; background-color: #7FB800; }
    .mmN { color: white; background-color: Salmon; }
    .deletion { color: #0D2C54; background-color: #0D2C54; }
    .deletion:active { display: none; }
    .insertion { color: cyan; font-weight: bold; text-decoration: underline; }
    .insertion:hover:after { content: attr(seq); background: #0D2C54; }
    .padding { opacity: 0; }
    .counts { color: #FF8C00; }
    .freq { color: #FF00FF; }
</style>
<script src="https://code.jquery.com/jquery-2.2.4.min.js" 
        integrity="sha256-BbhdlvQf/xTY9gja0Dq3HiwQF8LaCRTXxZKRutelT44=" 
        crossorigin="anonymous">
</script>
</head>
<body>
<div class="container">
    <div class="row">
        <button type="button" class="btn btn-info" id="expand_insertion">Expand All Insertions</button>
        <button type="button" class="btn btn-danger" id="remove_deletions">Remove All Deletions</button>
    </div>
</div>
<div class="alignment">
<table>
    <thead>
        <tr>
            <th>Alignment</th>
            <th>Counts</th>
            <th>Fraction</th>
        </tr>
    </thead>
    <tbody>
'''

HTML_TAIL = '''
        </tbody>
    </table>
</div>
<script>
 
$('#expand_insertion').click(function() {
 $('.insertion').each(function(idx, item) {
   let a = $(item);
   if(a.attr("seq") !== undefined) {
     a.html(a.html() + `<span style="color: blueviolet">${a.attr("seq")}</span>`);
     a.attr("seq", "");
   }
 });
$('#expand_insertion').prop("disabled", true);
});
 
$('#remove_deletions').click(function(){
 $('.deletion').remove();
 $('#remove_deletions').prop("disabled", true);
});
 
$('#remove_non_indel').click(function(){
 if($('#remove_non_indel').html().trim().substring(0,4) === "Hide") {
   $('.no_indel').hide(600);
   $('#remove_non_indel').html("Show Alignments with NO indels")
 } else {
   $('.no_indel').show(600);
   $('#remove_non_indel').html("Hide Alignments with NO indels")
 }
});
 
$('#remove_indel').click(function(){
 if($('#remove_indel').html().trim().substring(0,4) === "Hide") {
   $('.has_indel').hide(600);
   $('#remove_indel').html("Show Alignments with indels");
 } else {
   $('.has_indel').show(600);
   $('#remove_indel').html("Hide Alignments with indels");
 }
});
</script>
</body>
</html>
'''


def format_type(tag, frag):
    buffer = ""
    if tag == "M":
        buffer = ALIGN_MATCH.format(seq=frag)
    elif tag == "D":
        buffer = ALIGN_DELETION.format(seq='+' * len(frag))
    elif tag == "I":
        buffer = ALIGN_INSERTION.format(len=len(frag) - 1, seq=frag[1:], pre=frag[0])
    elif tag == "MM":
        for n in frag:
            buffer += ALIGN_MISMATCH.format(nt=n)
    return buffer


def df_to_html(df_alleles, ref, highlight, outfh, top_n=100):
    # print the template
    if highlight:
        ref_name, qw_name, qw, flank_bp = highlight.split(":")
        start, end = qw.split("-")
        outfh.write(ALIGN_TEMPLATE.format(prefix=ref[0:int(start) - 1], seq_to_highlight=ref[int(start) - 1:int(end)],
                                          postfix=ref[int(end):]))
    else:
        outfh.write(ALIGN_TEMPLATE.format(prefix=ref, seq_to_highlight="", postfix=""))

    # prepare the alignment counter
    aligncounter = Counter()
    total = 0
    # iterate through the reads
    for idx, row in df_alleles.iterrows():
        total += row["#Reads"]
        buffer = ""

        frag = ""
        tag = "M"
        for idx_c, c in enumerate(row["Aligned_Sequence"]):
            # a deletion
            if row["ref_positions"][idx_c] in row["all_deletion_positions"]:
                if tag != "D":
                    buffer += format_type(tag, frag)
                    frag = ""
                tag = "D"

            # a substitution
            elif row["ref_positions"][idx_c] in row["all_substitution_positions"]:
                if tag != "MM":
                    buffer += format_type(tag, frag)
                    frag = ""
                tag = "MM"

            # a insertion
            elif row["ref_positions"][idx_c] < 0:
                if row["ref_positions"][idx_c] == -1:
                    continue
                if tag != "I":
                    buffer += format_type(tag, frag[:-1])
                    frag = row["Aligned_Sequence"][idx_c - 1]
                tag = "I"

            # a nucleotide match
            else:
                if tag != "M":
                    buffer += format_type(tag, frag)
                    frag = ""
                tag = "M"
            frag += c

        if tag != "I":
            buffer += format_type(tag, frag)
        else:
            buffer += format_type("M", frag[0])
        aligncounter[buffer] += row["#Reads"]  # increment the counter

    # print out the top n
    for alignment, count in aligncounter.most_common(top_n):
        outfh.write(ALIGN_RECORD.format(alignment=alignment, count=count, frac=100 * count / total))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert CS2 Allele table to html")
    parser.add_argument("-f", "--crispresso_output_folder", dest="crispresso_output_folder", required=True,
                        help="crispresso_output_folder", type=str)
    parser.add_argument("-r", "--ref", dest="reference", required=True, help="reference name", type=str)
    parser.add_argument("-b", "--hl", dest="highlight", default="", help="reference position highlight", type=str)
    parser.add_argument("-n", "--topn", dest="topn", default=100, help="print the top N alignments", type=int)
    args = parser.parse_args()

    crispresso2_info = CRISPRessoShared.load_crispresso_info(args.crispresso_output_folder)
    z = zipfile.ZipFile(os.path.join(args.crispresso_output_folder,
                                     crispresso2_info['running_info']['allele_frequency_table_zip_filename']))
    zf = z.open(crispresso2_info['running_info']['allele_frequency_table_filename'])
    df_alleles = pd.read_csv(zf, sep="\t")
    df_alleles["all_deletion_positions"] = df_alleles["all_deletion_positions"].apply(eval)
    df_alleles["all_substitution_positions"] = df_alleles["all_substitution_positions"].apply(eval)
    df_alleles["ref_positions"] = df_alleles["ref_positions"].apply(eval)

    #    html_fh=open(os.path.join(args.crispresso_output_folder,args.reference+".html"),'w')
    html_fh = open(crispresso2_info['running_info']["name"] + "." + args.reference + ".html", 'w')
    html_fh.write(HTML_HEADER)
    df_to_html(df_alleles[df_alleles['Aligned_Reference_Names'] == args.reference],
               crispresso2_info["results"]["refs"][args.reference]["sequence"], args.highlight, html_fh)
    html_fh.write(HTML_TAIL)
    html_fh.close()
