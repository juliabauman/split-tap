#!/usr/bin/env bash
set -euo pipefail

IN_R1="$1"
IN_R2="$2"
OUT_R1="$3"
OUT_R2="$4"
MIN_READS="$5"
WARN_READS="$6"

# Count reads (gzipped FASTQ)
LINES_R1=$(zcat "$IN_R1" | wc -l | tr -d ' ')
if [[ "$LINES_R1" -eq 0 ]]; then
  echo "ERROR: $IN_R1 has 0 lines." >&2
  exit 1
fi

READS=$(( LINES_R1 / 4 ))

mkdir -p "$(dirname "$OUT_R1")"

# Emit warning if extremely small
if [[ "$READS" -lt "$WARN_READS" ]]; then
  echo "WARNING: CROP FASTQ extremely small ($READS reads < $WARN_READS)"
fi

# Passthrough if already large enough.
# Use a hardlink (not symlink): the input may be a Snakemake temp() file that
# gets deleted right after this job, which would leave a symlink dangling.
# A hardlink shares the inode, so the data survives; fall back to copy if the
# hardlink fails (e.g. cross-filesystem).
if [[ "$READS" -ge "$MIN_READS" ]]; then
  rm -f "$OUT_R1" "$OUT_R2"
  ln -f "$IN_R1" "$OUT_R1" 2>/dev/null || cp -f "$IN_R1" "$OUT_R1"
  ln -f "$IN_R2" "$OUT_R2" 2>/dev/null || cp -f "$IN_R2" "$OUT_R2"
  echo "boost_fastqs: passthrough ($READS >= $MIN_READS)"
  exit 0
fi

# Compute multiplier
MULT=$(( (MIN_READS + READS - 1) / READS ))
if [[ "$MULT" -lt 2 ]]; then
  MULT=2
fi

echo "boost_fastqs: boosting $READS → >=$MIN_READS (x$MULT)"

{
  for i in $(seq 1 "$MULT"); do
    zcat "$IN_R1"
  done
} | gzip -c > "$OUT_R1"

{
  for i in $(seq 1 "$MULT"); do
    zcat "$IN_R2"
  done
} | gzip -c > "$OUT_R2"
