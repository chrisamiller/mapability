#!/usr/bin/env bash
# Simulate, align and evaluate one chunk. Streams end to end (no FASTQ/BAM on disk).
#
# usage: run_chunk.sh -r ref.fa -c chunks.tsv -o OUTDIR [-i CHUNK_ID|INDEX] [-m normal|bisulfite]
#                     [-t THREADS] [-L 150] [-F 500] [-S] [-T 5] [-q 20] [-x "extra bwa-mem3 args"]
#
#   -L  read length; -F fragment (insert) length, outer end to outer end (paired only)
#   -S  single-end reads (one read per position; -F ignored)
#   -i  chunk id (e.g. chunk00012) or 1-based index into the distinct chunk ids.
#       Defaults to $LSB_JOBINDEX / $SLURM_ARRAY_TASK_ID (+1 if 0-based) if set.
#   -T  placement tolerance in bp (unclipped start vs. true origin)
#   -q  MAPQ threshold, only used for the per-chunk stats line
#
# Writes OUTDIR/<mode>_<layout>/chunks/<chunk_id>.gz and <chunk_id>.stats.tsv,
# where <layout> is pe<L>_f<F> or se<L> (so different settings never collide).
# A finished chunk is skipped on rerun, so array jobs can simply be resubmitted.
set -euo pipefail
# Find the bundled tools even when the runtime replaces PATH (e.g. LSF docker passes the host's PATH)
export PATH="/usr/local/bin:/opt/conda/bin:$PATH"

ref="" chunks="" out="" id="" mode="normal" threads=$(nproc) L=150 F=500 tol=5 minq=20 extra="" se=0
while getopts "r:c:o:i:m:t:L:F:ST:q:x:" opt; do
  case $opt in
    r) ref=$OPTARG ;; c) chunks=$OPTARG ;; o) out=$OPTARG ;; i) id=$OPTARG ;; m) mode=$OPTARG ;;
    t) threads=$OPTARG ;; L) L=$OPTARG ;; F) F=$OPTARG ;; T) tol=$OPTARG ;; q) minq=$OPTARG ;;
    x) extra=$OPTARG ;; S) se=1 ;;
    *) sed -n '2,16p' "$0"; exit 1 ;;
  esac
done
[[ -n $ref && -n $chunks && -n $out ]] || { sed -n '2,16p' "$0"; exit 1; }
[[ $mode == normal || $mode == bisulfite ]] || { echo "bad -m $mode" >&2; exit 1; }

if [[ -z $id ]]; then
  if [[ -n ${LSB_JOBINDEX:-} && ${LSB_JOBINDEX} != 0 ]]; then id=$LSB_JOBINDEX
  elif [[ -n ${SLURM_ARRAY_TASK_ID:-} ]]; then id=$((SLURM_ARRAY_TASK_ID + ${ARRAY_INDEX_OFFSET:-0}))
  else echo "no chunk given (-i) and no array index in environment" >&2; exit 1; fi
fi
if [[ $id =~ ^[0-9]+$ ]]; then
  id=$(cut -f1 "$chunks" | uniq | sed -n "${id}p")
  [[ -n $id ]] || { echo "chunk index out of range" >&2; exit 1; }
fi
grep -q "^${id}"$'\t' "$chunks" || { echo "chunk $id not in $chunks" >&2; exit 1; }

if (( se )); then layout=single tag=se$L; else layout=paired tag=pe${L}_f$F; fi
(( se || F >= L )) || { echo "fragment length (-F $F) must be >= read length (-L $L)" >&2; exit 1; }
dir=$out/${mode}_$tag/chunks
mkdir -p "$dir"
result=$dir/$id.gz
if [[ -s $result && -s $dir/$id.stats.tsv ]]; then
  echo "[run_chunk] $id already done, skipping" >&2; exit 0
fi

bwa_args=(-t "$threads" -v 1)
(( se )) || bwa_args+=(-p -I "$F,$((F / 10))")
if [[ $mode == bisulfite ]]; then
  bwa_args+=(--meth --meth-tags none)
fi
# shellcheck disable=SC2206
[[ -n $extra ]] && bwa_args+=($extra)

echo "[run_chunk] $(date '+%F %T') $id mode=$mode layout=$tag threads=$threads regions:" >&2
awk -v id="$id" '$1==id' "$chunks" | sed 's/^/    /' >&2

common=(--ref "$ref" --regions "$chunks" --chunk "$id" --read-len "$L" --frag-len "$F" --layout "$layout")
SECONDS=0
mapsim simulate "${common[@]}" --mode "$mode" \
  | bwa-mem3 mem "${bwa_args[@]}" "$ref" - \
  | mapsim evaluate "${common[@]}" --mode "$mode" --tolerance "$tol" --min-mapq "$minq" \
      -o "$result" --stats "$dir/$id.stats.tsv.tmp"
mv "$dir/$id.stats.tsv.tmp" "$dir/$id.stats.tsv"
echo "[run_chunk] $(date '+%F %T') $id done in ${SECONDS}s" >&2
cat "$dir/$id.stats.tsv" >&2
