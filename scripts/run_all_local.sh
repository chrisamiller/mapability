#!/usr/bin/env bash
# Run the whole pipeline on one machine (handy for small genomes / testing).
#
# usage: run_all_local.sh -r ref.fa -o OUTDIR [-m normal|bisulfite|both] [-s CHUNK_BP] [-t THREADS] [-q 20]
#                         [-L 150] [-F 500] [-S]
#   -L read length, -F fragment length (paired), -S single-end
#   ref.fa must already be prepared and indexed (prep_reference.sh).
set -euo pipefail
# Find the bundled tools even when the runtime replaces PATH (e.g. LSF docker passes the host's PATH)
export PATH="/usr/local/bin:/opt/conda/bin:$PATH"

ref="" out="" modes="both" size=50000000 threads=$(nproc) minq=20
L=150 F=500 se=""
while getopts "r:o:m:s:t:q:L:F:S" opt; do
  case $opt in
    r) ref=$OPTARG ;; o) out=$OPTARG ;; m) modes=$OPTARG ;; s) size=$OPTARG ;; t) threads=$OPTARG ;;
    q) minq=$OPTARG ;; L) L=$OPTARG ;; F) F=$OPTARG ;; S) se=-S ;; *) sed -n '2,7p' "$0"; exit 1 ;;
  esac
done
[[ -n $ref && -n $out ]] || { sed -n '2,7p' "$0"; exit 1; }
layout=(-L "$L" -F "$F" ${se:+"$se"})
[[ $modes == both ]] && modes="normal bisulfite"

mkdir -p "$out"
make_chunks.sh "$ref.fai" "$size" > "$out/chunks.tsv"
n=$(cut -f1 "$out/chunks.tsv" | uniq | wc -l)
for mode in $modes; do
  for i in $(seq 1 "$n"); do
    run_chunk.sh -r "$ref" -c "$out/chunks.tsv" -o "$out" -i "$i" -m "$mode" -t "$threads" -q "$minq" "${layout[@]}"
  done
  merge_chunks.sh -r "$ref" -c "$out/chunks.tsv" -o "$out" -m "$mode" -q "$minq" "${layout[@]}"
done
