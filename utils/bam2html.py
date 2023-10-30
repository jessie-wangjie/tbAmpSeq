"""
Smith-Waterman Aligner Wrapper

This module provide functions to traverse input FastA/Q file and align each
read to a template sequence with Smith-Waterman aligner, which can be found
in the clib module.
"""
import argparse
import os
from pysam import AlignmentFile
from collections import Counter

ALIGN_TEMPLATE = "<tr><td>" \
                 "<span class=\"template\">{prefix}</span>" \
                 "<span class=\"highlight\">{seq_to_highlight1}</span>" \
                 "<span class=\"template\">{middle}</span>" \
                 "<span class=\"highlight\">{seq_to_highlight2}</span>" \
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


def bam_to_html(bam_file, ref, fragment, highlight, outfh, top_n=100):
    # print the template
    h = []
    if fragment:
        ref_start, ref_end = fragment.split("-")
    else:
        ref_start = 1
        ref_end = len(ref)

    for k, window in enumerate(highlight):
        ref_name, qw_name, qw, flank_bp = window.split(":")
        start, end = qw.split("-")
        h.append((int(start), int(end)))
    h.sort()

    # ref = get_reference_sequence

    if len(h) == 1:
        outfh.write(ALIGN_TEMPLATE.format(prefix=ref[int(ref_start) - 1:h[0][0] - 1], seq_to_highlight1=ref[h[0][0] - 1:h[0][1]],
                                          middle="", seq_to_highlight2="", postfix=ref[h[0][1]:int(ref_end)]))
    if len(h) == 2:
        outfh.write(ALIGN_TEMPLATE.format(prefix=ref[int(ref_start) - 1:h[0][0] - 1], seq_to_highlight1=ref[h[0][0] - 1:h[0][1]],
                                          middle=ref[h[0][1]:h[1][0] - 1],
                                          seq_to_highlight2=ref[h[1][0] - 1:h[1][1]], postfix=ref[h[1][1]:int(ref_end)]))

    # prepare the alignment counter
    aligncounter = Counter()
    total = 0
    # iterate through the reads
    with AlignmentFile(bam_file, "rb") as bam:
        for alignment in bam:
            total += 1
            buffer = ""
            frag = ""

            if alignment.is_unmapped:
                continue
            # skip the bases are not in the fragment
            # if abs(row["ref_positions"][idx_c]) < int(ref_start) - 1 or abs(row["ref_positions"][idx_c]) > int(ref_end) - 1:
            #    continue
            if alignment.reference_start > 0:
                buffer += ALIGN_PADDING.format(seq='+' * alignment.reference_start)

            for qidx, ridx, base in alignment.get_aligned_pairs(with_seq=True):
                # a deletion
                if qidx is None:
                    buffer += ALIGN_DELETION.format(seq='+')
                # an insertion
                elif ridx is None:
                    frag += base
                    buffer += ALIGN_INSERTION.format(len=oplen, seq=base, pre=prev)
                # a mismatch
                elif base.islower():
                    buffer += ALIGN_MISMATCH.format(nt=base.upper())
                # a match
                else:
                    buffer += ALIGN_MATCH.format(seq=base)

            aligncounter[buffer] += 1  # increment the counter

    # print out the top n
    for alignment, count in aligncounter.most_common(top_n):
        outfh.write(ALIGN_RECORD.format(alignment=alignment, count=count, frac=100 * count / total))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert bam to html")
    parser.add_argument("-s", "--bam_file", dest="bam", required=True, help="bam file", type=str)
    parser.add_argument("-c", "--chr", dest="chr", required=True, help="chromosome to output", type=str)
    parser.add_argument("-o", "--output", dest="output", default=True, help="output file", type=str)
    parser.add_argument("-b", "--hl", dest="highlight", default=[], help="reference position highlight", action="append")
    parser.add_argument("-n", "--topn", dest="topn", default=10000000, help="print the top N alignments", type=int)
    parser.add_argument("-r", "--frag", dest="fragment", default="", help="specific region to plot", type=str)

    args = parser.parse_args()

    html_fh = open(args.output, 'w')
    html_fh.write(HTML_HEADER)
    bam_to_html(args.bam, args.chr, args.fragment, args.highlight, html_fh, args.topn)
    html_fh.write(HTML_TAIL)
    html_fh.close()
