#!/usr/bin/env python3
"""Build a small test genome with known mappability features.

usage: make_test_genome.py SOURCE.fa(.gz, bgzipped+faidx'd) OUT.fa

  chrA  = chr21:30,000,001-31,000,000 of the source, with a 1 kb N gap at 600,000
  chrB  = chr22:30,000,001-31,000,000 of the source, with
          - an exact copy of chrA[200,000:203,000) placed at chrB[500,000:503,000)  (segdup)
          - a copy of chrA[400,000:403,000) with every C->T placed at chrB[700,000:703,000)
            (distinct in normal space, identical to chrA after C->T conversion)
  chrUn_test1 = 2 kb unique sequence from chr20

Expected truth is written to OUT.truth.bed.
"""
import subprocess
import sys

src, out = sys.argv[1], sys.argv[2]


def fetch(region):
    seq = subprocess.run(["samtools", "faidx", src, region], check=True, capture_output=True, text=True).stdout
    return "".join(seq.split("\n")[1:]).upper()


a = list(fetch("chr21:30000001-31000000"))
b = list(fetch("chr22:30000001-31000000"))
u = fetch("chr20:40000001-40002000")

a[600000:601000] = "N" * 1000
b[500000:503000] = a[200000:203000]
b[700000:703000] = "".join(a[400000:403000]).replace("C", "T")

with open(out, "w") as fh:
    for name, s in (("chrA", "".join(a)), ("chrB", "".join(b)), ("chrUn_test1", u)):
        fh.write(f">{name}\n")
        for i in range(0, len(s), 60):
            fh.write(s[i:i + 60] + "\n")

with open(out + ".truth.bed", "w") as fh:
    fh.write("chrA\t200000\t203000\tsegdup_both_modes\n")
    fh.write("chrB\t500000\t503000\tsegdup_both_modes\n")
    fh.write("chrA\t600000\t601000\tN_gap\n")
    fh.write("chrA\t400000\t403000\tC2T_copy_bisulfite_OT_only\n")
    fh.write("chrB\t700000\t703000\tC2T_copy_bisulfite_OT_only\n")
