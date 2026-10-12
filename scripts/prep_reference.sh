#!/usr/bin/env bash
# Extract the contigs to score from a (b)gzipped or plain genome FASTA and
# build bwa-mem3 indexes.
#
# usage: prep_reference.sh -i genome.fa[.gz] -o OUTDIR [-c main|main+unplaced|all] [-l CONTIGS.txt]
#                          [-m none|normal|meth] [-t THREADS]
#
# Works for any species with UCSC (chr1, chrM) or Ensembl/GENCODE (1, MT) style names:
#   -c main            numbered chromosomes, X, Y, M/MT (with or without "chr")
#      main+unplaced   (default) plus unplaced/unlocalized scaffolds: chr*_random, chrUn_*,
#                      and accession-named scaffolds like GL456210.1 / KI270728.1
#                      (alts, decoys, HLA, EBV and patches are excluded)
#      all             every contig in the input
#   -l FILE            keep exactly the contigs listed in FILE (one name per line); overrides -c.
#                      Use for other naming schemes (e.g. RefSeq NC_000067.7) or W/Z/2L-style names.
#   -m normal          build the standard index only (~76 GB RAM for hg38)
#      meth            (default) build standard + bisulfite seed index (~152 GB RAM for hg38);
#                      this one index serves both modes
#      none            only write ref.fa (index later)
set -euo pipefail
# Find the bundled tools even when the runtime replaces PATH (e.g. LSF docker passes the host's PATH)
export PATH="/usr/local/bin:/opt/conda/bin:$PATH"

in="" out="" contigs="main+unplaced" list="" index="meth" threads=$(nproc)
while getopts "i:o:c:l:m:t:" opt; do
  case $opt in
    i) in=$OPTARG ;; o) out=$OPTARG ;; c) contigs=$OPTARG ;; l) list=$OPTARG ;; m) index=$OPTARG ;;
    t) threads=$OPTARG ;;
    *) sed -n '2,19p' "$0"; exit 1 ;;
  esac
done
[[ -n $in && -n $out ]] || { sed -n '2,19p' "$0"; exit 1; }
mkdir -p "$out"
ref=$out/ref.fa

main_pat='^(chr)?([0-9]+|X|Y|M|MT)$'
unplaced_pat='^chr[0-9XY]+_[A-Za-z0-9]+_random$|^chrUn_[A-Za-z0-9]+$|^[A-Z]{2}[0-9]{6}\.[0-9]+$'
case $contigs in
  main)          pat=$main_pat ;;
  main+unplaced) pat="$main_pat|$unplaced_pat" ;;
  all)           pat='.' ;;
  *) echo "unknown -c $contigs" >&2; exit 1 ;;
esac

if [[ ! -s $ref.fai ]]; then
  # list contig names (works for plain or gzipped input)
  names=$(mktemp)
  if [[ -s $in.fai ]]; then cut -f1 "$in.fai"; else
    { [[ $in == *.gz ]] && gzip -dc "$in" || cat "$in"; } | grep '^>' | sed 's/^>//; s/[[:space:]].*//'
  fi > "$names.all"
  if [[ -n $list ]]; then
    wanted=$(tr -d '\r' < "$list" | grep -v '^[[:space:]]*$' || true)
    missing=$(grep -vxFf "$names.all" <<< "$wanted" || true)
    [[ -z $missing ]] || { echo "[prep] contigs in $list not found in input:" $missing >&2; exit 1; }
    grep -xFf <(echo "$wanted") "$names.all" > "$names" || true
  else
    grep -E "$pat" "$names.all" > "$names" || true
  fi
  [[ -s $names ]] || { echo "no contigs selected (pattern: ${list:-$pat})" >&2; exit 1; }
  echo "[prep] keeping $(wc -l < "$names") of $(wc -l < "$names.all") contigs, e.g.:" \
    $(head -5 "$names") >&2
  dropped=$(grep -vxFf "$names" "$names.all" || true)
  [[ -z $dropped ]] || echo "[prep] dropping $(wc -l <<< "$dropped") contigs, e.g.:" $(head -5 <<< "$dropped") >&2
  rm -f "$names.all"

  if [[ -s $in.fai ]] && { [[ $in != *.gz ]] || [[ -s $in.gzi ]]; }; then
    samtools faidx -r "$names" "$in" > "$ref.tmp"
  else
    # plain gzip (not bgzip): stream and filter
    { [[ $in == *.gz ]] && gzip -dc "$in" || cat "$in"; } |
      awk -v list="$names" 'BEGIN{while((getline n < list)>0) keep[n]=1}
           /^>/{name=substr($1,2); p=(name in keep)} p' > "$ref.tmp"
  fi
  mv "$ref.tmp" "$ref"
  samtools faidx "$ref"
  rm -f "$names"
fi

case $index in
  none) ;;
  normal) [[ -s $ref.bwt.2bit.64 ]] || bwa-mem3 index -t "$threads" "$ref" ;;
  meth)   [[ -s $ref.meth.bwt.2bit.64 || -n $(ls "$ref".meth.*bwt* 2>/dev/null) ]] || bwa-mem3 index --meth -t "$threads" "$ref" ;;
  *) echo "unknown -m $index" >&2; exit 1 ;;
esac
awk '{n++; s+=$2} END{printf "[prep] %s: %d contigs, %d bp\n", FILENAME, n, s}' "$ref.fai" >&2
