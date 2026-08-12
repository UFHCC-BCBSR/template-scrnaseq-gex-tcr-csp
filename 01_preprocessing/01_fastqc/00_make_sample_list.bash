#!/bin/bash
#
# Build sample_fastqs.tsv by finding every R1/R2 fastq pair in the delivery
# directory named in config/config.env.
#
# Sequencing cores nest fastqs in deep, inconsistent folder trees, so this
# searches the whole tree rather than assuming a layout. Files are matched by
# their library tag (GEX, CSP, TCR, BCR) as it appears in the filename.
#
# Run:
#     bash 00_make_sample_list.bash
#
# Then check the output before submitting the FastQC array job:
#     wc -l sample_fastqs.tsv
#     head sample_fastqs.tsv

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/load_config.bash"
check_config

OUT_TSV="${SCRIPT_DIR}/sample_fastqs.tsv"

if [ ! -d "$RAW_FASTQ_DIR" ]; then
  echo "ERROR: RAW_FASTQ_DIR does not exist: $RAW_FASTQ_DIR" >&2
  exit 1
fi

echo "Searching for fastqs under: $RAW_FASTQ_DIR"
echo "Library types: $LIBRARY_TYPES"
echo

: > "$OUT_TSV"

total=0
for tag in $LIBRARY_TYPES; do
  n=0
  # Find R1 files for this library tag, then derive the matching R2.
  while IFS= read -r -d '' r1; do
    r2="${r1/_R1_/_R2_}"
    if [ ! -f "$r2" ]; then
      echo "  WARNING: no R2 mate for $(basename "$r1"), skipping" >&2
      continue
    fi
    # Label = the fastq prefix, e.g. "1-GEX" from 1-GEX_S1_L001_R1_001.fastq.gz
    label="$(basename "$r1" | sed 's/_S[0-9]*_L[0-9]*_R1_001\.fastq\.gz$//')"
    printf '%s\t%s\t%s\n' "$label" "$r1" "$r2" >> "$OUT_TSV"
    n=$((n + 1))
  done < <(find -L "$RAW_FASTQ_DIR" -type f -name "*${tag}*_R1_001.fastq.gz" -print0 | sort -z)

  echo "  ${tag}: ${n} read pairs"
  total=$((total + n))
done

# Keep the file deterministic regardless of filesystem ordering.
sort -o "$OUT_TSV" "$OUT_TSV"

echo
echo "Wrote ${total} read pairs to ${OUT_TSV}"

if [ "$total" -eq 0 ]; then
  echo
  echo "ERROR: no fastqs matched. Check that:" >&2
  echo "  - RAW_FASTQ_DIR points at the right directory" >&2
  echo "  - your filenames contain the library tags in LIBRARY_TYPES" >&2
  echo "  - your files end in _R1_001.fastq.gz / _R2_001.fastq.gz" >&2
  exit 1
fi

echo
echo "Compare these counts against your sequencing core's delivery manifest."
echo "A library sequenced across N lanes should show N read pairs per pool."
echo "Then submit FastQC with:  sbatch 01_fastqc.sbatch"
