#!/bin/bash
#
# Run the per-sample R pipeline for one analysis group, or for every group.
#
#     bash run_all.bash              # every group in config/samples.tsv
#     bash run_all.bash TISSUE_A     # one group
#
# Steps run in order and each depends on the previous one, so the script stops
# at the first failure rather than carrying on with stale inputs.
#
# Steps 07 (publish) and 09 (Loupe) are optional and are allowed to fail
# without stopping the run, since they depend on external configuration and
# an extra R package respectively.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/load_config.bash"

cd "${SCRIPT_DIR}/R"

if [ $# -gt 0 ]; then
  GROUPS_TO_RUN="$*"
else
  GROUPS_TO_RUN="$(sample_groups)"
fi

REQUIRED_STEPS=(
  01_load_data.R
  02_qc_filtering.R
  03_vdj_integration.R
  04_norm_and_cluster.R
  05_cell_annotation.R
  06_de.R
)

OPTIONAL_STEPS=(
  07_prepare_report.R
  09_loupe.R
)

for group in $GROUPS_TO_RUN; do
  echo "############################################################"
  echo "# Group: ${group}"
  echo "############################################################"

  for step in "${REQUIRED_STEPS[@]}"; do
    echo
    echo "--- ${step} (${group}) ---"
    Rscript "$step" "$group"
  done

  for step in "${OPTIONAL_STEPS[@]}"; do
    echo
    echo "--- ${step} (${group}, optional) ---"
    Rscript "$step" "$group" || echo "NOTE: ${step} did not complete; continuing."
  done

  echo
  echo "--- 08_report.Rmd (${group}) ---"
  Rscript -e "rmarkdown::render('08_report.Rmd', params = list(dataset = '${group}'), output_file = paste0('report_', '${group}', '.html'))" \
    || echo "NOTE: report rendering did not complete; continuing."
done

echo
echo "All requested groups finished."
echo "Next, combine groups with:  cd ../integrated_analysis/R && Rscript 01_integrate_data.R"
