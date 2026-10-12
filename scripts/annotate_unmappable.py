#!/usr/bin/env python3
"""Annotate unmappable regions with the likely reason (repeat class, segdup, N gap).

usage:
  annotate_unmappable.py --bed unmappable.bed --nbed N.bed \
      --rmsk rmsk.txt.gz [--simple simpleRepeat.txt.gz] [--segdup genomicSuperDups.txt.gz] \
      [--alias chromAlias.txt.gz] [--genome ref.fa.fai] [--min-frac 0.5] [--list-min-frac 0.1] --out PREFIX

Inputs are UCSC database tables (hgdownload.soe.ucsc.edu/goldenPath/<db>/database/). Contig
names are translated to the BED's naming (e.g. chr1 -> 1, chr1_GL456210_random -> GL456210.1)
using chromAlias when given.

Each region gets the fraction of its bases covered by each category. The primary reason is the
highest-priority category covering >= --min-frac of the region:

  N > segdup > satellite > simple_repeat > low_complexity > LINE > SINE > LTR > DNA > other_repeat

otherwise "unexplained". Repeats are reported per family (e.g. LINE/L1, SINE/B1(Alu)).
The annotations column lists every label covering >= --list-min-frac of the region, in priority
order ("." if none), e.g. "segdup,LTR/ERVK".

Outputs:
  PREFIX.annotated.bed   chrom start end primary_reason length annotations   (0-based, half-open)
  PREFIX.summary.tsv     counts of regions and bp per primary reason, plus bp-level overlap and
                         genome-wide background for each category
"""
import argparse
import bisect
import collections
import gzip
import sys

PRIORITY = ["N", "segdup", "satellite", "simple_repeat", "low_complexity",
            "LINE", "SINE", "LTR", "DNA", "other_repeat"]

# SINE family display names (mouse B1 is labeled family "Alu" by UCSC RepeatMasker)
SINE_FAMILY = {"Alu": "B1(Alu)", "B2": "B2", "B4": "B4", "ID": "ID", "MIR": "MIR"}


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def log(msg):
    print(f"[annotate] {msg}", file=sys.stderr, flush=True)


def load_alias(path, target_names):
    """Map UCSC contig names -> names used in the BED."""
    if not path:
        return lambda chrom: chrom if chrom in target_names else None
    m = {}
    for line in opener(path):
        alias, ucsc = line.rstrip("\n").split("\t")[:2]
        if alias in target_names:
            m[ucsc] = alias
    for name in target_names:  # names already in UCSC form
        m.setdefault(name, name)
    return m.get


def rmsk_category(rclass, rfamily):
    """Return (category, label) for a RepeatMasker class/family."""
    cls = rclass.rstrip("?")
    if cls == "Satellite":
        return "satellite", "satellite/" + rfamily
    if cls == "Simple_repeat":
        return "simple_repeat", "simple_repeat"
    if cls == "Low_complexity":
        return "low_complexity", "low_complexity"
    if cls == "LINE":
        return "LINE", "LINE/" + rfamily.rstrip("?")
    if cls == "SINE":
        fam = rfamily.rstrip("?")
        return "SINE", "SINE/" + SINE_FAMILY.get(fam, "other")
    if cls == "LTR":
        return "LTR", "LTR/" + rfamily.rstrip("?")
    if cls == "DNA":
        return "DNA", "DNA/" + rfamily.rstrip("?")
    return "other_repeat", "other_repeat/" + rclass


class Track:
    """Per-contig sorted, merged intervals for one label."""

    def __init__(self):
        self.raw = collections.defaultdict(list)
        self.starts = {}
        self.ends = {}

    def add(self, chrom, start, end):
        self.raw[chrom].append((start, end))

    def finalize(self):
        for chrom, interval in self.raw.items():
            interval.sort()
            ms, me = [], []
            for start, end in interval:
                if ms and start <= me[-1]:
                    if end > me[-1]:
                        me[-1] = end
                else:
                    ms.append(start)
                    me.append(end)
            self.starts[chrom], self.ends[chrom] = ms, me
        self.raw = None

    def overlap(self, chrom, start, end):
        st = self.starts.get(chrom)
        if not st:
            return 0
        en = self.ends[chrom]
        i = bisect.bisect_right(en, start)  # first interval ending after start
        t = 0
        while i < len(st) and st[i] < end:
            t += min(end, en[i]) - max(start, st[i])
            i += 1
        return t

    def total(self):
        return sum(end - start for chrom in self.starts
                   for start, end in zip(self.starts[chrom], self.ends[chrom]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True, help="regions to annotate (e.g. *.unmappable.bed)")
    ap.add_argument("--nbed", help="reference N runs (e.g. *.N.bed)")
    ap.add_argument("--rmsk", required=True, help="UCSC rmsk.txt.gz")
    ap.add_argument("--simple", help="UCSC simpleRepeat.txt.gz (TRF); merged into simple_repeat")
    ap.add_argument("--segdup", help="UCSC genomicSuperDups.txt.gz")
    ap.add_argument("--alias", help="UCSC chromAlias.txt.gz, to translate contig names")
    ap.add_argument("--genome", help="contig lengths (name<TAB>length, e.g. .fai) for genome-wide background")
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="fraction of a region a category must cover to be the primary reason (0.5)")
    ap.add_argument("--list-min-frac", type=float, default=0.1,
                    help="fraction of a region a label must cover to appear in the annotations column (0.1)")
    ap.add_argument("--out", required=True, help="output prefix")
    a = ap.parse_args()

    regions = []
    for line in opener(a.bed):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        regions.append((f[0], int(f[1]), int(f[2])))
    names = {r[0] for r in regions}
    genome = {}
    if a.genome:
        for line in opener(a.genome):
            f = line.split("\t")
            genome[f[0]] = int(f[1])
        names |= set(genome)
    log(f"{len(regions):,} regions on {len(names)} contigs")
    to_target = load_alias(a.alias, names)

    tracks = collections.defaultdict(Track)  # label -> Track
    label_cat = {}
    unmapped_contigs = collections.Counter()

    def add(label, cat, chrom, start, end):
        t = to_target(chrom)
        if t is None:
            unmapped_contigs[chrom] += 1
            return
        tracks[label].add(t, start, end)
        label_cat[label] = cat

    if a.nbed:
        for line in opener(a.nbed):
            f = line.split("\t")
            add("N", "N", f[0], int(f[1]), int(f[2]))
    log("reading RepeatMasker")
    for line in opener(a.rmsk):
        f = line.split("\t")
        cat, label = rmsk_category(f[11], f[12])
        add(label, cat, f[5], int(f[6]), int(f[7]))
    if a.simple:
        log("reading TRF simple repeats")
        for line in opener(a.simple):
            f = line.split("\t")
            add("simple_repeat", "simple_repeat", f[1], int(f[2]), int(f[3]))
    if a.segdup:
        log("reading segmental duplications")
        for line in opener(a.segdup):
            f = line.split("\t")
            add("segdup", "segdup", f[1], int(f[2]), int(f[3]))
    if unmapped_contigs:
        log(f"skipped annotations on {len(unmapped_contigs)} contigs not in the BED/genome "
            f"(e.g. {', '.join(list(unmapped_contigs)[:3])})")

    # category-level tracks = union of their labels
    cat_track = collections.defaultdict(Track)
    for label, t in tracks.items():
        for chrom, interval in t.raw.items():
            cat_track[label_cat[label]].raw[chrom].extend(interval)
    any_rep = Track()
    for cat in PRIORITY[2:]:
        if cat in cat_track:
            for chrom, interval in cat_track[cat].raw.items():
                any_rep.raw[chrom].extend(interval)
    for t in list(tracks.values()) + list(cat_track.values()) + [any_rep]:
        t.finalize()
    labels = sorted(tracks, key=lambda l: (PRIORITY.index(label_cat[l]), l))
    log(f"{len(labels)} annotation labels")

    # annotate
    by_reason_n = collections.Counter()
    by_reason_bp = collections.Counter()
    cat_bp = collections.Counter()
    sine_bp = collections.Counter()
    any_rep_bp = 0
    total_bp = 0
    with open(a.out + ".annotated.bed", "w") as out:
        out.write("#chrom\tstart\tend\tprimary_reason\tlength\tannotations\n")
        for chrom, start, end in regions:
            n = end - start
            total_bp += n
            cf = {cat: cat_track[cat].overlap(chrom, start, end) for cat in cat_track}
            for cat, bp in cf.items():
                cat_bp[cat] += bp
            any_rep_bp += any_rep.overlap(chrom, start, end)
            primary = next((cat for cat in PRIORITY if cf.get(cat, 0) >= a.min_frac * n), "unexplained")
            parts = []
            for l in labels:
                bp = tracks[l].overlap(chrom, start, end)
                if bp:
                    parts.append((l, bp))
                    if label_cat[l] == "SINE":
                        sine_bp[l] += bp
            # primary label = most-covering label within the primary category
            if primary != "unexplained":
                primary = max((p for p in parts if label_cat[p[0]] == primary), key=lambda p: p[1])[0]
            ann = ",".join(l for l, bp in parts if bp >= a.list_min_frac * n) or "."
            out.write(f"{chrom}\t{start}\t{end}\t{primary}\t{n}\t{ann}\n")
            key = primary.split("/")[0] if not primary.startswith("SINE/") else primary
            by_reason_n[key] += 1
            by_reason_bp[key] += n

    # summary
    gtot = sum(genome.values()) if genome else None
    order = (["N", "segdup", "satellite", "simple_repeat", "low_complexity", "LINE"]
             + sorted(k for k in by_reason_n if k.startswith("SINE/"))
             + ["LTR", "DNA", "other_repeat", "unexplained"])
    with open(a.out + ".summary.tsv", "w") as fh:
        fh.write("# A. primary reason per region (each region counted once; priority order, "
                 f">= {a.min_frac:.0%} of region)\n")
        fh.write("primary_reason\tregions\tpct_regions\tbp\tpct_bp\n")
        nreg = len(regions)
        for k in order:
            if by_reason_n[k]:
                fh.write(f"{k}\t{by_reason_n[k]}\t{100 * by_reason_n[k] / nreg:.2f}\t{by_reason_bp[k]}\t"
                         f"{100 * by_reason_bp[k] / total_bp:.2f}\n")
        fh.write(f"total\t{nreg}\t100.00\t{total_bp}\t100.00\n\n")

        fh.write("# B. bp-level overlap: fraction of all unmappable bp covered by each category "
                 "(categories overlap, so this does not sum to 100%)")
        fh.write("; genome_pct = same category's share of the whole genome; enrichment = ratio\n" if gtot else "\n")
        fh.write("category\tunmappable_bp\tpct_of_unmappable" + ("\tgenome_bp\tpct_of_genome\tenrichment\n" if gtot else "\n"))
        rows = [(cat, cat_bp[cat], cat_track[cat]) for cat in PRIORITY if cat in cat_track]
        rows += [(l, sine_bp[l], tracks[l]) for l in labels if label_cat[l] == "SINE"]
        rows += [("any_repeat", any_rep_bp, any_rep)]
        for name, bp, tr in rows:
            fh.write(f"{name}\t{bp}\t{100 * bp / total_bp:.2f}")
            if gtot:
                g = tr.total()
                pu, pg = bp / total_bp, g / gtot
                fh.write(f"\t{g}\t{100 * pg:.2f}\t{pu / pg if pg else float('nan'):.2f}")
            fh.write("\n")
    log(f"wrote {a.out}.annotated.bed and {a.out}.summary.tsv")


if __name__ == "__main__":
    main()
