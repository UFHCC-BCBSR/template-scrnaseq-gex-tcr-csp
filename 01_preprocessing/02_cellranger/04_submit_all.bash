#!/bin/bash
#
# Submit one cellranger multi job per pool.
#
# Each job is submitted from inside its own <POOL>_run directory, because
# cellranger multi writes its output tree into the working directory.
#
# Run:
#     bash 04_submit_all.bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/load_config.bash"
check_config

POOLS="$(sample_pools)"
submitted=0

for POOL in $POOLS; do
  run="${SCRIPT_DIR}/${POOL}_run"
  job="${run}/cellranger_${POOL}.sbatch"

  if [ ! -f "$job" ]; then
    echo "ERROR: $job not found. Run 'bash 03_make_per_pool_runs.bash' first." >&2
    exit 1
  fi

  # Skip pools that already have completed output, so re-running this script
  # after a partial failure does not clobber good results.
  if [ -f "${run}/${POOL}_multi/outs/per_sample_outs" ] || [ -d "${run}/${POOL}_multi/outs" ]; then
    echo "  ${POOL}: output already exists at ${run}/${POOL}_multi/outs, skipping"
    echo "           (delete that directory to re-run this pool)"
    continue
  fi

  echo -n "  ${POOL}: "
  ( cd "$run" && sbatch "cellranger_${POOL}.sbatch" )
  submitted=$((submitted + 1))
done

echo
echo "Submitted ${submitted} job(s)."
echo
echo "Watch them with:      squeue -u \$USER"
echo "Follow a log with:    tail -f <POOL>_run/logs/cellranger_<POOL>.*.out"
echo "When done, check:     bash 05_check_runs.bash"
