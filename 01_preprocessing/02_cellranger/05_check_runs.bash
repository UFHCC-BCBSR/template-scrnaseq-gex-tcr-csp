#!/bin/bash
#
# Check whether each cellranger multi run finished, and point at the outputs
# you should look at.
#
# Run:
#     bash 05_check_runs.bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/load_config.bash"

POOLS="$(sample_pools)"
n_ok=0
n_bad=0

for POOL in $POOLS; do
  run="${SCRIPT_DIR}/${POOL}_run"
  outs="${run}/${POOL}_multi/outs"
  expected=$(samples_data | awk -F'\t' -v p="$POOL" '$1==p' | wc -l | tr -d ' ')

  echo "=== ${POOL} ==="

  if [ ! -d "$outs" ]; then
    echo "  NOT FINISHED - no outs directory yet."
    if compgen -G "${run}/logs/*.err" > /dev/null; then
      echo "  Last lines of the most recent error log:"
      tail -n 15 "$(ls -t "${run}"/logs/*.err | head -1)" | sed 's/^/    /'
    fi
    n_bad=$((n_bad + 1))
    echo
    continue
  fi

  ws="${outs}/per_sample_outs"
  if [ -d "$ws" ]; then
    got=$(find "$ws" -maxdepth 1 -mindepth 1 -type d | wc -l | tr -d ' ')
    echo "  Samples recovered: ${got} (expected ${expected})"
    if [ "$got" -ne "$expected" ]; then
      echo "  WARNING: sample count does not match config/samples.tsv." >&2
      echo "           Check the hashtag assignments and the web summary." >&2
    fi
  fi

  echo "  Web summaries to open in a browser:"
  find "$outs" -name web_summary.html | sed 's/^/    /'
  n_ok=$((n_ok + 1))
  echo
done

echo "-----------------------------------------------"
echo "Finished: ${n_ok}    Not finished / failed: ${n_bad}"
echo
echo "Open each web_summary.html and check:"
echo "  - estimated number of cells is plausible, not near zero"
echo "  - the multiplexing section shows the expected samples per pool"
echo "  - no large red warning banners"
