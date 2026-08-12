#!/bin/bash
#
# Run FastQC on every fastq pair listed in sample_fastqs.tsv, as a SLURM array
# job (one array task per pair).
#
# Generate the input list first:
#     bash 00_make_sample_list.bash
#
# Then submit. The array size is read from the file, so submit via sbatch with
# an explicit array range:
#     sbatch --array=1-$(wc -l < sample_fastqs.tsv) 01_fastqc.sbatch
#
# Or use the wrapper, which works this out for you:
#     bash 02_submit_fastqc.bash
#
#SBATCH --job-name=fastqc
#SBATCH --output=logs/fastqc_%A_%a.out
#SBATCH --error=logs/fastqc_%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
source "${SCRIPT_DIR}/../../config/load_config.bash"

SAMPLE_LIST="${SCRIPT_DIR}/sample_fastqs.tsv"
OUTDIR="${SCRIPT_DIR}/reports"
mkdir -p "$OUTDIR" "${SCRIPT_DIR}/logs"

if [ ! -f "$SAMPLE_LIST" ]; then
  echo "ERROR: $SAMPLE_LIST not found. Run 00_make_sample_list.bash first." >&2
  exit 1
fi

if [ "${USE_MODULES:-true}" = "true" ]; then
  module load "${FASTQC_MODULE:-fastqc}"
fi

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
if [ -z "$LINE" ]; then
  echo "ERROR: no line ${SLURM_ARRAY_TASK_ID} in $SAMPLE_LIST" >&2
  exit 1
fi

SAMPLE=$(echo "$LINE" | cut -f1)
FQ1=$(echo "$LINE" | cut -f2)
FQ2=$(echo "$LINE" | cut -f3)

echo "Running FastQC for $SAMPLE:"
echo "  R1: $FQ1"
echo "  R2: $FQ2"

fastqc -t "${SLURM_CPUS_PER_TASK:-2}" -o "$OUTDIR" "$FQ1" "$FQ2"

echo "Done. Reports in $OUTDIR"
