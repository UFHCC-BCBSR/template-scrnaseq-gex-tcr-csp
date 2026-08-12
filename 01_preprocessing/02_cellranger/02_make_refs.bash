#!/bin/bash
#
# Download the Cell Ranger references for the species set in config/config.env.
#
# This only needs to be done once per species, ever. The downloads are large
# (roughly 10-30 GB depending on species) and can take a while.
#
# If your cluster already provides shared 10x references, you do not need this
# script at all: point REF_GEX / REF_VDJ in 03_make_per_pool_runs.bash at the
# shared copies instead, or symlink them into refs/.
#
# Run:
#     bash 02_make_refs.bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/load_config.bash"
check_config

REFS_DIR="${SCRIPT_DIR}/refs"
mkdir -p "$REFS_DIR"
cd "$REFS_DIR"

# ---- reference URLs per species ---------------------------------------------
# From https://www.10xgenomics.com/support/software/cell-ranger/downloads
# Check that page for newer builds; update these URLs as needed.

case "$SPECIES" in
  mouse)
    GEX_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
    GEX_NAME="refdata-gex-GRCm39-2024-A"
    VDJ_URL="https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0.tar.gz"
    VDJ_NAME="refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0"
    ;;
  human)
    GEX_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
    GEX_NAME="refdata-gex-GRCh38-2024-A"
    VDJ_URL="https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0.tar.gz"
    VDJ_NAME="refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0"
    ;;
  *)
    echo "ERROR: unsupported SPECIES '${SPECIES}' (expected mouse or human)" >&2
    echo "       For another species, add its reference URLs to this script." >&2
    exit 1
    ;;
esac

# Does this experiment need a V(D)J reference at all?
NEED_VDJ=false
if has_library TCR || has_library BCR; then
  NEED_VDJ=true
fi

echo "Species: ${SPECIES}"
echo "Refs directory: ${REFS_DIR}"
echo

# ---- 1. gene expression reference -------------------------------------------
if [ -d "$GEX_NAME" ]; then
  echo "GEX reference already present: ${GEX_NAME} (skipping)"
else
  echo "Downloading GEX reference: ${GEX_NAME}"
  wget --continue "$GEX_URL"
  tar -xzf "$(basename "$GEX_URL")"
  rm -f "$(basename "$GEX_URL")"
fi

# ---- 2. V(D)J reference ------------------------------------------------------
if [ "$NEED_VDJ" = false ]; then
  echo "No TCR or BCR library configured; skipping V(D)J reference."
elif [ -d "$VDJ_NAME" ]; then
  echo "VDJ reference already present: ${VDJ_NAME} (skipping)"
else
  echo "Downloading VDJ reference: ${VDJ_NAME}"
  wget --continue "$VDJ_URL"
  tar -xzf "$(basename "$VDJ_URL")"
  rm -f "$(basename "$VDJ_URL")"
fi

echo
echo "References ready:"
echo "  GEX: ${REFS_DIR}/${GEX_NAME}"
if [ "$NEED_VDJ" = true ]; then
  echo "  VDJ: ${REFS_DIR}/${VDJ_NAME}"
fi
