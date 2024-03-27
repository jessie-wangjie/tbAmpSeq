"""
Smith-Waterman Aligner Wrapper

This module provide functions to traverse input FastA/Q file and align each
read to a template sequence with Smith-Waterman aligner, which can be found
in the clib module.
"""
import argparse
import re
from pysam import AlignmentFile, FastaFile
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
ALIGN_RECORD = "<tr class=\"{indel}\"><td>{alignment}</td><td>{count}</td><td>{frac:.2f}%</td></tr>"

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
        <button type="button" class="btn btn-primary" id="remove_non_indel">Hide Alignments with NO indels</button>
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


def bam_to_html(bam_file, fasta, region, highlight, outfh, top_n):
    # open reference file
    fa = FastaFile(fasta)

    # parse region
    if ":" in region:
        m = re.match(r"(.*):(\d+)-(\d+)", region)
        chr = m.group(1)
        ref_start = int(m.group(2))
        ref_end = int(m.group(3))
    else:
        chr = region
        ref_start = 1
        ref_end = fa.get_reference_length(chr)
    ref_seq = fa.fetch(chr, ref_start - 1, ref_end)

    # print the template
    h = []
    for qw in highlight:
        start, end = qw.split("-")
        h.append((int(start) - ref_start + 1, int(end) - ref_start + 1))
    h.sort()

    if len(h) == 1:
        outfh.write(ALIGN_TEMPLATE.format(prefix=ref_seq[ref_start - 1:h[0][0] - 1], seq_to_highlight1=ref_seq[h[0][0] - 1:h[0][1]],
                                          middle="", seq_to_highlight2="", postfix=ref_seq[h[0][1]:ref_end]))
    elif len(h) == 2:
        outfh.write(ALIGN_TEMPLATE.format(prefix=ref_seq[ref_start - 1:h[0][0] - 1], seq_to_highlight1=ref_seq[h[0][0] - 1:h[0][1]],
                                          middle=ref_seq[h[0][1]:h[1][0] - 1],
                                          seq_to_highlight2=ref_seq[h[1][0] - 1:h[1][1]], postfix=ref_seq[h[1][1]:ref_end]))
    else:
        outfh.write(ALIGN_TEMPLATE.format(prefix=ref_seq, seq_to_highlight1="", middle="", seq_to_highlight2="", postfix=""))

    # prepare the alignment counter
    aligncounter = Counter()
    total = 0

    # iterate through the reads
    bam = AlignmentFile(bam_file, "rb")
    # only fetch the alignments overlapping with the shown region
    for alignment in bam.fetch(chr, int(ref_start) - 1, int(ref_end)):
        total += 1
        buffer = ""

        # skip the unmapped reads or not primary alignments or mapQ < 1
        if alignment.is_unmapped or alignment.is_secondary or alignment.mapping_quality < 1:
            continue

        # add paddings for alignments which don't start 1bp of shown region
        if alignment.reference_start + 1 > int(ref_start):
            buffer += ALIGN_PADDING.format(seq='+' * (alignment.reference_start - int(ref_start) + 1))

        i = alignment.query_alignment_start
        pos = alignment.get_aligned_pairs(with_seq=True)
        indel = "no_indel"

        while i < len(pos):
            # skip bases which are not within the shown region
            if (pos[i][1] is not None) and (pos[i][1] < int(ref_start) - 1 or pos[i][1] >= int(ref_end)):
                i = i + 1
                continue

            # 3' soft-clip
            if (pos[i][1] is None) and (pos[i][2] is None) and (int(alignment.query_alignment_end) <= pos[i][0]):
                i = i + 1
                continue

            # a deletion
            if pos[i][0] is None:
                if h[0][0] - 1 <= pos[i][1] <= h[0][1]:
                    indel = "has_indel"
                buffer += ALIGN_DELETION.format(seq='+')
            # an insertion
            elif pos[i][1] is None:
                prev_pos = pos[i - 1][1]
                if h[0][0] - 1 <= prev_pos <= h[0][1]:
                    indel = "has_indel"

                insertions = alignment.query_sequence[pos[i][0]]
                while i < len(pos) - 1 and pos[i + 1][1] is None:
                    i = i + 1
                    insertions += alignment.query_sequence[pos[i][0]]
                if int(ref_start) - 1 <= prev_pos < int(ref_end):
                    prev = buffer[-1]
                    if prev == ">":
                        prev = ""
                    else:
                        buffer = buffer[:-1]
                    buffer += ALIGN_INSERTION.format(len=len(insertions), seq=insertions, pre=prev)
            # mismatch
            elif pos[i][2].islower():
                buffer += ALIGN_MISMATCH.format(nt=alignment.query_sequence[pos[i][0]].upper())
            # match
            else:
                buffer += ALIGN_MATCH.format(seq=pos[i][2])
            i = i + 1
        aligncounter[tuple([buffer, indel])] += 1  # increment the counter

    # print out the top n
    for alignment, count in aligncounter.most_common(top_n):
        outfh.write(ALIGN_RECORD.format(alignment=alignment[0], count=count, frac=100 * count / total, indel=alignment[1]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert bam to html")
    parser.add_argument("-s", "--bam_file", dest="bam", required=True, help="bam file", type=str)
    parser.add_argument("-f", "--reference_file", dest="reference", required=True, help="reference fasta file", type=str)
    parser.add_argument("-r", "--region", dest="region", required=True, help="specific region to plot", type=str)
    parser.add_argument("-o", "--output", dest="output", default=True, help="output file", type=str)
    parser.add_argument("-b", "--hl", dest="highlight", default=[], help="reference position highlight", action="append")
    parser.add_argument("-n", "--topn", dest="topn", default=10000000, help="print the top N alignments", type=int)

    args = parser.parse_args()

    html_fh = open(args.output, 'w')
    html_fh.write(HTML_HEADER)
    bam_to_html(args.bam, args.reference, args.region, args.highlight, html_fh, args.topn)
    html_fh.write(HTML_TAIL)
    html_fh.close()
