#!/bin/bash
#
# Combine all the individual FastQC reports into one MultiQC report.
#
# Run after the FastQC array job has finished:
#     bash 03_multiqc.bash
#
# Then open reports/multiqc_report.html in a browser. See the README section
# "What to check, and what is a red flag" for how to read it.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/load_config.bash"

REPORT_DIR="${SCRIPT_DIR}/reports"

if [ ! -d "$REPORT_DIR" ]; then
  echo "ERROR: $REPORT_DIR not found. Run the FastQC job first." >&2
  exit 1
fi

n_zip=$(find "$REPORT_DIR" -maxdepth 1 -name '*_fastqc.zip' | wc -l | tr -d ' ')
if [ "$n_zip" -eq 0 ]; then
  echo "ERROR: no *_fastqc.zip files in $REPORT_DIR." >&2
  echo "       Did the array job finish? Check logs/ for errors." >&2
  exit 1
fi

echo "Summarising ${n_zip} FastQC reports..."

if [ "${USE_MODULES:-true}" = "true" ]; then
  # MultiQC is commonly available as its own module; adjust if yours differs.
  module load multiqc 2>/dev/null || echo "NOTE: no multiqc module; assuming it is on PATH"
fi

multiqc "$REPORT_DIR" --outdir "$REPORT_DIR" --force

echo
echo "Wrote ${REPORT_DIR}/multiqc_report.html"
echo "Open it in a browser and review it before running Cell Ranger."
