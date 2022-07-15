#!/usr/bin/env python
"""
On-Target analysis for Amp-seq and 3Primer-Seq
Information from Benchling

functions need by tbOnT pipeline
"""

import subprocess
import re


def align_primer(seq, index, chromosome, adapter=""):
    seq = seq.upper()
    adapter = adapter.upper()

    if adapter in seq:
        seq = seq.replace(adapter, "")

    fastq = ">seq" + "\n" + seq + "\n+" + "\n" + seq
    bwa_out = subprocess.Popen("echo -e '%s' | bwa fastmap %s -" % (fastq, index), stdout=subprocess.PIPE, shell=True)
    pos = subprocess.check_output(["grep", "EM"], stdin=bwa_out.stdout).split()[4:]
    for p in pos:
        m = re.match(r"(.*):([+|-])(\d+)", p.decode())
        if chromosome == m.group(1):
            return {"chr": m.group(1), "start": int(m.group(3)), "end": int(m.group(3)) + len(seq) - 1,
                    "strand": m.group(2), "seq": seq}


def get_cut_site(seq, guide):
    # return 1-index
    # cut is always the left of the cutting site
    guide = guide.replace("U", "T")
    if guide in seq:
        p5 = seq.find(guide) + 1
        p3 = p5 + len(guide) - 1
        cut = p3 - 3
        guide_strand = "+"
    else:
        p3 = seq.find(reverse_complement(guide)) + 1
        p5 = p3 + len(guide) - 1
        cut = p3 - 1 + 3
        guide_strand = "-"
    return {"5P": p5, "3P": p3, "cut": cut, "strand": guide_strand, "seq": guide}


def get_seq(twobit_file, chromosome, start, end, strand):
    seq = subprocess.check_output(
        "twoBitToFa -seq=%s -start=%s -end=%s %s stdout | grep -v \> | xargs | sed 's/ //g'" % (
            chromosome, start - 1, end, twobit_file), shell=True).decode().rstrip()
    if strand == "-":
        seq = reverse_complement(seq)
    return seq.upper()


def get_beacon_seq(seq1, sp1_strand, seq2="", sp2_strand=""):
    # beacon1 and beacon2 should have at least 5nt overlap.
    beacon = seq1.upper()
    if sp1_strand == "+":
        beacon = reverse_complement(seq1)
    if seq2 != "":
        beacon2 = seq2.upper()
        if sp2_strand == "+":
            beacon2 = reverse_complement(seq2)
        idx = beacon.find(beacon2[0:5])
        beacon = beacon[0:idx] + beacon2
    return beacon


def reverse_complement(seq):
    nt_complement = dict({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '_': '_', '-': '-', 'U': 'A'})
    return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])
