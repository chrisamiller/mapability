#!/usr/bin/env bash
# Check all chunks finished, then merge them into the final outputs.
#
# usage: merge_chunks.sh -r ref.fa -c chunks.tsv -o OUTDIR [-m normal|bisulfite] [-L 150] [-F 500] [-S] [-q 20]
#   -L/-F/-S must match what run_chunk.sh was given (they select OUTDIR/<mode>_<layout>/)
#
# Outputs in OUTDIR/<mode>_<layout>/  (<layout> = pe<L>_f<F> or se<L>):
#   <name>.unmappable.bed        merged intervals with best MAPQ < q (0-based, half-open);
#                                N bases count as mappable if a good read covers them
#   <name>.bestmapq.bedgraph.gz  best MAPQ of any correctly placed read; -1 = none
#   <name>.single_strand.bed     bisulfite only: bases covered by good reads from just one of
#                                OT/OB; 4th column = the strand that maps (OT_only / OB_only)
#   <name>.N.bed                 reference N runs
#   <name>.stats.tsv             per-chunk read placement stats
# where <name> = mapability.<mode>.<layout>.q<q> (unless -n is given)
set -euo pipefail
# Find the bundled tools even when the runtime replaces PATH (e.g. LSF docker passes the host's PATH)
export PATH="/usr/local/bin:/opt/conda/bin:$PATH"

ref="" chunks="" out="" mode="normal" minq=20 name="" L=150 F=500 se=0
while getopts "r:c:o:m:q:n:L:F:S" opt; do
  case $opt in
    r) ref=$OPTARG ;; c) chunks=$OPTARG ;; o) out=$OPTARG ;; m) mode=$OPTARG ;; q) minq=$OPTARG ;;
    n) name=$OPTARG ;; L) L=$OPTARG ;; F) F=$OPTARG ;; S) se=1 ;;
    *) sed -n '2,15p' "$0"; exit 1 ;;
  esac
done
[[ -n $ref && -n $chunks && -n $out ]] || { sed -n '2,15p' "$0"; exit 1; }
if (( se )); then tag=se$L; else tag=pe${L}_f$F; fi
name=${name:-mapability.$mode.$tag.q$minq}
dir=$out/${mode}_$tag
list=$dir/chunk_files.txt

missing=0
: > "$list"
for id in $(cut -f1 "$chunks" | uniq); do
  if [[ -s $dir/chunks/$id.gz && -s $dir/chunks/$id.stats.tsv ]]; then
    echo "$dir/chunks/$id.gz" >> "$list"
  else
    echo "missing: $id" >&2; missing=$((missing + 1))
  fi
done
(( missing == 0 )) || { echo "$missing chunk(s) not finished; not merging" >&2; exit 1; }

mapsim merge --ref "$ref" --regions "$list" --out "$dir/$name" --min-mapq "$minq"

{ head -1 "$(head -1 "$list" | sed 's/\.gz$/.stats.tsv/')"
  sed 's/\.gz$/.stats.tsv/' "$list" | xargs -n 200 tail -q -n +2; } > "$dir/$name.stats.tsv"
awk -F'\t' 'NR>1{p+=$3; u+=$4; w+=$5; c+=$6; g+=$7}
  END{printf "[merge] primary=%d unmapped=%d misplaced=%d correct=%d correct_mapq_ok=%d (%.3f%%)\n",
      p,u,w,c,g, p?100*g/p:0}' "$dir/$name.stats.tsv" >&2
echo "[merge] wrote $dir/$name.unmappable.bed" >&2
