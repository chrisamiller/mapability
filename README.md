# mapability

Exhaustive paired-end or single-end mappability of a genome with
[bwa-mem3](https://github.com/fg-labs/bwa-mem3), in normal and bisulfite modes.

For every fragment start position, a perfect 2×150 bp pair with a 500 bp insert (defaults) is simulated, aligned
back with bwa-mem3, and checked against its known origin. A base is **unmappable** if no correctly placed
read with MAPQ ≥ 20 covers it. Reference N bases are treated like any other base: an N covered by a
good read (e.g. at the edge of a gap, aligned as a mismatch) is mappable; gap interiors that no read
can reach are unmappable.

## Method

* **Contigs** (`prep_reference.sh -c`): numbered chromosomes, X, Y, M/MT plus unplaced/unlocalized scaffolds
  (no alts, decoys, HLA, EBV, patches). Works with UCSC names (`chr1`, `chrM`, `chr1_GL456210_random`,
  `chrUn_JH584304`) and Ensembl/GENCODE names (`1`, `MT`, `JH584295.1`, `KI270728.1`), so human, mouse and
  other species work unchanged. `-c main` drops the scaffolds, `-c all` keeps everything, and
  `-l contigs.txt` keeps exactly the listed names (for other schemes, e.g. RefSeq `NC_000067.7`).
* **Reads**: read names encode origin (`contig:fragstart0:strand`). Pairs where *both* reads are entirely N are
  skipped (they cannot align); reads that are partly N are kept.
  * normal: R1 = ref[f, f+150), R2 = revcomp(ref[f+350, f+500)). Only one strand is simulated: the
    minus-strand pair is the same pair with R1/R2 swapped.
  * bisulfite (directional library, every C converted): OT pair (R1 = C→T of the top strand) and OB pair
    (R1 = C→T of the bottom strand), aligned with `bwa-mem3 mem --meth` (collapsed scoring, bwameth-like).
* **Read layout** (`-L`, `-F`, `-S` on `run_chunk.sh`, `merge_chunks.sh`, `run_all_local.sh`; `READ_LEN`,
  `FRAG_LEN`, `SINGLE=1` for `lsf/submit_lsf.sh`):
  * `-L` read length (default 150).
  * `-F` fragment/insert length, outer end of R1 to outer end of R2 (default 500, i.e. a 200 bp unsequenced
    gap between the reads; TLEN = 500 in the SAM).
  * `-S` single-end: one read of length L per position (normal), or one OT and one OB read (bisulfite).
* **Alignment**: paired: `bwa-mem3 mem -p -I F,F/10` (fixed insert distribution so results do not depend on
  chunking); single-end: `bwa-mem3 mem`.
* **Correct placement**: primary record, origin contig and strand, unclipped start within 5 bp of origin.
  Each correctly placed read credits the bases it aligns to with its MAPQ; a base's score is the max.
  Because each read is scored on its own MAPQ, a read in a repeat whose mate anchors in unique sequence
  counts as mappable (bwa pairing gives it a high MAPQ) — i.e. this is paired-end mappability. Use `-S`
  for single-end mappability, which flags considerably more sequence.
* Simulation → alignment → evaluation is streamed; no FASTQ/BAM is written.

## Outputs (`OUT/<mode>_<layout>/`, layout = `pe<L>_f<F>` or `se<L>`, e.g. `bisulfite_pe150_f500`)

| file | content |
|---|---|
| `mapability.<mode>.<layout>.q20.unmappable.bed` | merged unmappable intervals (best MAPQ < 20). Standard BED: 0-based, half-open, so 1-based positions 2,3,4 are `chr 1 4` |
| `...bestmapq.bedgraph.gz` | best MAPQ of any correctly placed read per base; `-1` = none. Re-threshold without re-aligning: `zcat x.bedgraph.gz \| awk '$4<30' \| bedtools merge` |
| `...single_strand.bed` | bisulfite only: bases covered by good reads from just one strand (OT or OB). 4th column is the strand that maps, `OT_only` / `OB_only`. These are mappable in the main BED, but methylation on the other strand's C's can't be measured there |
| `...N.bed` | reference N runs (for reference; not forced into the unmappable BED) |
| `...stats.tsv` | per-chunk counts: primary reads, unmapped, misplaced, correct, correct with MAPQ ≥ q |

## Running

```bash
docker buildx build --platform linux/amd64 -t yourrepo/mapability:1.4.0 --push .

# 1. reference + indexes (once; the --meth index serves both modes)
prep_reference.sh -i b38.fa.gz -o work/ref -c main+unplaced -m meth -t 32

# 2a. on a cluster: see lsf/submit_lsf.sh (array of run_chunk.sh jobs + dependent merge_chunks.sh)
# 2b. or one machine:
run_all_local.sh -r work/ref/ref.fa -o work/out -m both -t 32
run_all_local.sh -r work/ref/ref.fa -o work/out -m both -t 32 -L 100 -F 350   # other read layout
run_all_local.sh -r work/ref/ref.fa -o work/out -m both -t 32 -S              # single-end 150
```

Steps individually:

```bash
make_chunks.sh ref.fa.fai 50000000 > chunks.tsv          # small contigs are packed together
run_chunk.sh -r ref.fa -c chunks.tsv -o out -m normal -i 7 -t 16   # -i defaults to $LSB_JOBINDEX
merge_chunks.sh -r ref.fa -c chunks.tsv -o out -m normal -q 20     # give the same -L/-F/-S as run_chunk.sh
```

`run_chunk.sh` skips chunks that are already finished, so a failed array can just be resubmitted.
Pass extra aligner options with `-x "..."`.

## Annotating unmappable regions

`annotate_unmappable.py` labels each unmappable region with the likely reason, using UCSC
RepeatMasker, TRF simple repeats and segmental duplications. UCSC contig names are translated to
yours (e.g. `chr1` → `1`, `chr1_GL456210_random` → `GL456210.1`) with UCSC's chromAlias table.

```bash
fetch_ucsc_annotations.sh mm10 annotations/mm10        # GRCm38 = mm10, GRCm39 = mm39, GRCh38 = hg38
annotate_unmappable.py --bed out/normal_pe150_f500/mapability.normal.pe150_f500.q20.unmappable.bed \
    --nbed out/normal_pe150_f500/mapability.normal.pe150_f500.q20.N.bed \
    --rmsk annotations/mm10/rmsk.txt.gz --simple annotations/mm10/simpleRepeat.txt.gz \
    --segdup annotations/mm10/genomicSuperDups.txt.gz --alias annotations/mm10/chromAlias.txt.gz \
    --genome work/ref/ref.fa.fai --out out/normal_pe150_f500/mapability.normal.pe150_f500.q20
```

(In the container with a host `PATH`, call it as `/usr/local/bin/python3 /usr/local/bin/annotate_unmappable.py`.)

* `*.annotated.bed`: `chrom start end primary_reason length annotations`, e.g.
  `1  6431630  6434607  segdup  2977  segdup,LTR/ERVK`.
  * `primary_reason`: the highest-priority category covering ≥50% of the region (`--min-frac`), shown at
    family level. Priority: N > segdup > satellite > simple_repeat > low_complexity > LINE > SINE > LTR >
    DNA > other_repeat; otherwise `unexplained`, which is mostly regions where no single category reaches 50%.
  * `annotations`: every label covering ≥10% of the region (`--list-min-frac`), in priority order; `.` if none.
  * SINEs are split by family: `SINE/B1(Alu)` (mouse B1 is UCSC's "Alu" family), `B2`, `B4`, `ID`, `MIR`.
* `*.summary.tsv`: (A) regions and bp per primary reason; (B) bp-level overlap of each category with all
  unmappable bp, its share of the whole genome (`--genome`), and the enrichment ratio.

## Resources (hg38, up to 32 cores / 256 GB per job)

| step | cores | memory | notes |
|---|---|---|---|
| `prep_reference.sh -m meth` | 32 | 200 GB | bwa-mem3 docs: ~152 GB for the dual index (~76 GB for `-m normal`); ~35 GB disk |
| normal chunks (50 Mb) | 16 | 32 GB | index resident ~10.5 GB; ~62 jobs |
| bisulfite chunks (25 Mb) | 16 | 48 GB | seed+original index ~22 GB; twice the pairs per bp; ~125 jobs |
| merge | 2 | 16 GB | a byte per genome position per strand type (~3.1 GB normal, ~6.2 GB bisulfite) |

16-core jobs are recommended over 32-core ones: bwa-mem3 scales nearly linearly, so total CPU time is
about the same, but two jobs fit per 32-core node and schedule more easily, and the per-job index load is
small next to chunk runtime. Use 32 cores / 64 GB (normal) or 96 GB (bisulfite) if job slots are scarce.
The single-threaded `mapsim simulate`/`evaluate` stages run at ~0.75-0.9 M pairs/s, well ahead of
bwa-mem3 even at 32 threads.

Time: on a small test (4 Apple-silicon threads) normal mode ran ~65 k pairs/s and bisulfite ~25 k
(OT+OB) pairs/s. Whole-genome repeats will make real throughput lower; plan on roughly 100-300
CPU-hours for normal and 2-3× that for bisulfite, and calibrate by running one chunk first and reading
the `done in Ns` line from its log.

## Test

```bash
python3 test/make_test_genome.py ~/annotations/b38.fa.gz test/genome.fa   # needs samtools
docker run --rm -v $PWD/test:/data mapability:1.4.0 bash -c \
  'prep_reference.sh -i /data/genome.fa -o /data/ref -c all && run_all_local.sh -r /data/ref/ref.fa -o /data/out -s 400000 -t 4'
```

Expected (both modes): `chrA 200496 202504`, `chrB 500496 502504` (interior of a 3 kb exact segdup;
the outer ~500 bp is rescued by mates in unique sequence) and the N gap minus the few edge N's
that reads overhang (`chrA 600004 600996` normal, `chrA 600009 600991` bisulfite). A 3 kb copy
that differs only by C→T (`chrA:400000-403000` vs `chrB:700000-703000`) is *not* flagged in
bisulfite mode: OT reads are ambiguous there (MAPQ 0) but OB reads still place uniquely, so both
copies appear in `single_strand.bed` as `OB_only`.

With `-L 100 -F 300` the segdup interior is `chrA 200296 202704` (mates rescue ~300 bp); with `-S`
(single-end) it is `chrA 200140 202860`, and several short natural repeats are flagged that paired
reads resolve.
