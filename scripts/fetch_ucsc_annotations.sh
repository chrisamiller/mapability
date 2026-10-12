#!/usr/bin/env bash
# Download the UCSC tables used by annotate_unmappable.py.
#
# usage: fetch_ucsc_annotations.sh UCSC_DB OUTDIR      e.g. mm10, mm39, hg38
#
# Fetches rmsk (RepeatMasker), simpleRepeat (TRF), genomicSuperDups (segmental duplications)
# and chromAlias (to translate chr1 -> 1, chr1_GL456210_random -> GL456210.1, ...).
# Check that UCSC_DB matches the assembly you aligned to (e.g. GRCm38 = mm10, GRCm39 = mm39).
set -euo pipefail
export PATH="/usr/local/bin:/opt/conda/bin:$PATH"
db=${1:?usage: fetch_ucsc_annotations.sh UCSC_DB OUTDIR}
out=${2:?usage: fetch_ucsc_annotations.sh UCSC_DB OUTDIR}
mkdir -p "$out"
base=https://hgdownload.soe.ucsc.edu/goldenPath/$db/database
for t in rmsk simpleRepeat genomicSuperDups chromAlias; do
  if [[ -s $out/$t.txt.gz ]]; then echo "[fetch] $t present" >&2; continue; fi
  echo "[fetch] $base/$t.txt.gz" >&2
  if ! curl -sfL -o "$out/$t.txt.gz.tmp" "$base/$t.txt.gz"; then
    rm -f "$out/$t.txt.gz.tmp"
    if [[ $t == genomicSuperDups || $t == chromAlias ]]; then
      echo "[fetch] warning: $t not available for $db; continuing without it" >&2; continue
    fi
    echo "[fetch] error: could not download $t for $db" >&2; exit 1
  fi
  mv "$out/$t.txt.gz.tmp" "$out/$t.txt.gz"
done
date -u '+fetched %Y-%m-%dT%H:%M:%SZ from '"$base" > "$out/SOURCE.txt"
