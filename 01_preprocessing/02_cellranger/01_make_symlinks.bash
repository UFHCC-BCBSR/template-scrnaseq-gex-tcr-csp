#!/bin/bash
#
# Consolidate the sequencing delivery into one flat directory of symlinks per
# library type, so each cellranger multi run can point --fastqs at a single dir.
#
# Sequencing cores deliver fastqs in deep, nested, inconsistent folder trees.
# Cell Ranger wants a flat directory. Symlinks give us that without copying
# hundreds of gigabytes.
#
# Files are located by their library tag in the filename (e.g. 1-GEX_...,
# 2-CSP_...), which is robust to whatever folder layout the core used.
#
# Run:
#     bash 01_make_symlinks.bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/load_config.bash"
check_config

OUTDIR="${SCRIPT_DIR}/symlinks"

if [ ! -d "$RAW_FASTQ_DIR" ]; then
  echo "ERROR: RAW_FASTQ_DIR does not exist: $RAW_FASTQ_DIR" >&2
  exit 1
fi

# Map each library tag to the symlink subdirectory Cell Ranger will read.
# These names must match the [libraries] paths written by 03_make_per_pool_runs.bash.
libdir_for_tag() {
  case "$1" in
    GEX) echo "gex_fastqs"  ;;
    CSP) echo "csp_fastqs"  ;;
    TCR) echo "vdjt_fastqs" ;;
    BCR) echo "vdjb_fastqs" ;;
    *)   echo "" ;;
  esac
}

echo "Delivery root: $RAW_FASTQ_DIR"
echo "Linking into:  $OUTDIR"
echo

for tag in $LIBRARY_TYPES; do
  subdir="$(libdir_for_tag "$tag")"
  if [ -z "$subdir" ]; then
    echo "WARNING: unknown library type '${tag}', skipping" >&2
    continue
  fi

  dest="${OUTDIR}/${subdir}"
  mkdir -p "$dest"

  # Link every R1/R2 fastq for this library type, anywhere under the delivery.
  find -L "$RAW_FASTQ_DIR" -type f -name "*-${tag}_*_R[12]_001.fastq.gz" -print0 \
    | while IFS= read -r -d '' fq; do
        ln -sf "$fq" "$dest/"
      done

  n=$(find "$dest" -maxdepth 1 -type l -name '*.fastq.gz' | wc -l | tr -d ' ')
  echo "  ${tag}: linked ${n} fastq.gz -> ${dest}"

  if [ "$n" -eq 0 ]; then
    echo "    WARNING: nothing matched *-${tag}_*_R[12]_001.fastq.gz" >&2
    echo "    Check that your filenames use this tag. Adjust the find pattern" >&2
    echo "    above if your core uses a different naming convention." >&2
    echo "    Example filenames from the delivery:" >&2
    find -L "$RAW_FASTQ_DIR" -type f -name '*.fastq.gz' 2>/dev/null \
      | head -3 | sed 's|.*/|      |' >&2
  fi
done

echo
echo "Compare these counts against your delivery manifest before continuing."
echo "A library sequenced across N lanes gives 2N files (R1 + R2) per pool."
echo "If the numbers look wrong, stop and investigate rather than running Cell Ranger."
