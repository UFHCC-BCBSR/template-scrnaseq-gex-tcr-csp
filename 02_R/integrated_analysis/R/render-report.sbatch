#!/bin/bash
#
# Render an R Markdown report on a compute node.
#
# Usage:
#     sbatch render-report.sbatch <report.Rmd> [output.html]
#
# Examples:
#     sbatch render-report.sbatch 08_report.Rmd
#     sbatch render-report.sbatch 08_report.Rmd report_TISSUE_A.html
#
# Submit from the directory containing the .Rmd file. Account and QoS are read
# from config/config.env and passed on the command line, because SLURM parses
# #SBATCH directives before any shell code runs.
#
# If your allocation is not set up, submit manually instead:
#     sbatch --account=ACCT --qos=QOS render-report.sbatch 08_report.Rmd
#
#SBATCH --job-name=render_report
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50gb
#SBATCH --time=03:00:00
#SBATCH --output=logs/render_%j.log
#SBATCH --error=logs/render_%j.err

set -euo pipefail

REPORT="${1:-}"
OUTPUT="${2:-}"

if [ -z "$REPORT" ]; then
  echo "ERROR: no report specified." >&2
  echo "Usage: sbatch render-report.sbatch <report.Rmd> [output.html]" >&2
  exit 1
fi

if [ ! -f "$REPORT" ]; then
  echo "ERROR: report not found: $REPORT" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../config/load_config.bash"

mkdir -p logs

if [ "${USE_MODULES:-true}" = "true" ]; then
  module purge
  module load "${R_MODULE:-R/4.5}"
  module load "${PANDOC_MODULE:-pandoc}"
fi

# rmarkdown needs to know where pandoc lives.
if command -v pandoc > /dev/null 2>&1; then
  RSTUDIO_PANDOC="$(dirname "$(command -v pandoc)")"
  export RSTUDIO_PANDOC
  echo "Using pandoc: $(command -v pandoc)"
else
  echo "WARNING: pandoc not found on PATH; rendering may fail." >&2
fi

# R_LIBS_USER is the user's personal library; R sets it by default. Do not
# hard-code a path here, it differs per user and per R version.
echo "R library paths:"
Rscript -e '.libPaths()'

if [ -n "$OUTPUT" ]; then
  Rscript -e "rmarkdown::render('${REPORT}', output_file = '${OUTPUT}')"
else
  Rscript -e "rmarkdown::render('${REPORT}')"
fi

echo "Rendered ${REPORT}"
