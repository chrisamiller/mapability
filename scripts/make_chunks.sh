#!/usr/bin/env bash
# Split the reference into chunks of fragment start positions.
# Small contigs are packed together so no job is tiny.
#
# usage: make_chunks.sh ref.fa.fai CHUNK_BP > chunks.tsv
# output columns: chunk_id  contig  start0  end0   (one chunk may span several lines/contigs)
set -euo pipefail
# Find the bundled tools even when the runtime replaces PATH (e.g. LSF docker passes the host's PATH)
export PATH="/usr/local/bin:/opt/conda/bin:$PATH"
fai=${1:?usage: make_chunks.sh ref.fa.fai [CHUNK_BP]}
size=${2:-50000000}

awk -v size="$size" 'BEGIN{OFS="\t"; id=0; acc=0}
{
  len=$2
  for (s=0; s<len; s+=size) {
    e=s+size; if (e>len) e=len
    if (acc>0 && acc+(e-s) > size) { id++; acc=0 }
    print sprintf("chunk%05d", id), $1, s, e
    acc+=e-s
    if (acc>=size) { id++; acc=0 }
  }
}' "$fai"
