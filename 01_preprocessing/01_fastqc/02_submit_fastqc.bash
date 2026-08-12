#!/bin/bash
#
# Submit the FastQC array job with the right array size and SLURM account.
#
# SLURM reads #SBATCH directives before any shell code runs, so the account,
# QoS and array range cannot come from config.env inside the job script itself.
# This wrapper passes them on the sbatch command line instead.
#
# Run:
#     bash 02_submit_fastqc.bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/load_config.bash"
check_config

SAMPLE_LIST="${SCRIPT_DIR}/sample_fastqs.tsv"

if [ ! -f "$SAMPLE_LIST" ]; then
  echo "ERROR: $SAMPLE_LIST not found." >&2
  echo "       Run 'bash 00_make_sample_list.bash' first." >&2
  exit 1
fi

N=$(wc -l < "$SAMPLE_LIST" | tr -d ' ')
if [ "$N" -eq 0 ]; then
  echo "ERROR: $SAMPLE_LIST is empty." >&2
  exit 1
fi

mkdir -p "${SCRIPT_DIR}/logs" "${SCRIPT_DIR}/reports"

echo "Submitting FastQC for ${N} fastq pairs..."

sbatch \
  --account="$SLURM_ACCOUNT" \
  --qos="$SLURM_QOS" \
  --array="1-${N}" \
  --chdir="$SCRIPT_DIR" \
  "${SCRIPT_DIR}/01_fastqc.sbatch"

echo
echo "Watch progress with:  squeue -u \$USER"
echo "When all tasks finish, summarise them with:  bash 03_multiqc.bash"
