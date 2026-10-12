#!/usr/bin/env bash
# Example LSF submission (job arrays + dependent merge). Edit the settings block for your cluster.
#
# usage: submit_lsf.sh normal|bisulfite
#   read layout via env: READ_LEN (150), FRAG_LEN (500, paired only), SINGLE=1 for single-end
#
# Prerequisite (one job, run once; the --meth index serves both modes):
#   bsub -q $QUEUE -G $GROUP -a "docker($IMAGE)" -n 32 -R 'span[hosts=1] rusage[mem=200GB]' -M 200GB \
#        -o $WORK/logs/prep.%J.log \
#        /usr/local/bin/prep_reference.sh -i /path/to/genome.fa.gz -o $WORK/ref -c main+unplaced -m meth -t 32
set -euo pipefail
mode=${1:?usage: submit_lsf.sh normal|bisulfite}

# ---- settings -------------------------------------------------------------
IMAGE=${IMAGE:-yourrepo/mapability:1.4.0}
QUEUE=${QUEUE:-general}
GROUP=${GROUP:-compute-yourlab}
WORK=${WORK:-/scratch/yourlab/mapability}           # must be visible inside the container
export LSF_DOCKER_VOLUMES="${LSF_DOCKER_VOLUMES:-$WORK:$WORK}"
REF=$WORK/ref/ref.fa
MINQ=20
READ_LEN=${READ_LEN:-150}
FRAG_LEN=${FRAG_LEN:-500}
SINGLE=${SINGLE:-0}
layout=(-L "$READ_LEN" -F "$FRAG_LEN")
(( SINGLE )) && layout+=(-S)
MAX_RUNNING=40                                       # array throttle
if [[ $mode == normal ]]; then
  CHUNK_BP=50000000; CORES=16; MEM=32GB
else
  CHUNK_BP=25000000; CORES=16; MEM=48GB              # 2x pairs per bp, bigger seed index
fi
# ---------------------------------------------------------------------------

out=$WORK/out
mkdir -p "$out" "$WORK/logs"
chunks=$out/chunks.$mode.tsv

# make the chunk table (cheap; done here on the head node if mapsim scripts are on PATH,
# otherwise as a small blocking job)
if command -v make_chunks.sh >/dev/null; then
  make_chunks.sh "$REF.fai" "$CHUNK_BP" > "$chunks"
else
  bsub -K -q "$QUEUE" -G "$GROUP" -a "docker($IMAGE)" -o "$WORK/logs/chunks.%J.log" \
    bash -c "/usr/local/bin/make_chunks.sh $REF.fai $CHUNK_BP > $chunks"
fi
n=$(cut -f1 "$chunks" | uniq | wc -l)
echo "submitting $n $mode chunks"

if (( SINGLE )); then tag=se$READ_LEN; else tag=pe${READ_LEN}_f$FRAG_LEN; fi
name=map_${mode}_$tag
bsub -q "$QUEUE" -G "$GROUP" -a "docker($IMAGE)" \
  -J "${name}[1-$n]%$MAX_RUNNING" -n "$CORES" -R "span[hosts=1] rusage[mem=$MEM]" -M "$MEM" \
  -o "$WORK/logs/$name.%J.%I.log" \
  /usr/local/bin/run_chunk.sh -r "$REF" -c "$chunks" -o "$out" -m "$mode" -t "$CORES" -q "$MINQ" "${layout[@]}"

bsub -q "$QUEUE" -G "$GROUP" -a "docker($IMAGE)" \
  -J "${name}_merge" -w "done(${name})" -n 2 -R "rusage[mem=16GB]" -M 16GB \
  -o "$WORK/logs/${name}_merge.%J.log" \
  /usr/local/bin/merge_chunks.sh -r "$REF" -c "$chunks" -o "$out" -m "$mode" -q "$MINQ" "${layout[@]}"
