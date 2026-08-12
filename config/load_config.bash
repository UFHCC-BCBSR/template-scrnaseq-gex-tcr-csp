#!/bin/bash
#
# Shared config loader for the bash / sbatch scripts.
#
# Usage, from any script in the repo:
#     source "$(git rev-parse --show-toplevel)/config/load_config.bash"
# or, without relying on git:
#     source "/path/to/repo/config/load_config.bash"
#
# After sourcing you have every KEY from config/config.env as a shell variable,
# plus REPO_ROOT, CELLRANGER_DIR, SAMPLES_TSV and the helper functions below.

# Resolve the repo root from this file's own location, so it works no matter
# which directory the calling script was launched from.
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$CONFIG_DIR")"
CONFIG_FILE="${CONFIG_DIR}/config.env"
SAMPLES_TSV="${CONFIG_DIR}/samples.tsv"

CELLRANGER_DIR="${REPO_ROOT}/01_preprocessing/02_cellranger"
FASTQC_DIR="${REPO_ROOT}/01_preprocessing/01_fastqc"

if [ ! -f "$CONFIG_FILE" ]; then
  echo "ERROR: config file not found: $CONFIG_FILE" >&2
  return 1 2>/dev/null || exit 1
fi

# shellcheck disable=SC1090
set -a
source "$CONFIG_FILE"
set +a

# ---- validation -------------------------------------------------------------
# Fail early and loudly rather than producing empty paths deep inside a job.

config_die() {
  echo "ERROR: $*" >&2
  echo "       Fix it in ${CONFIG_FILE}" >&2
  return 1 2>/dev/null || exit 1
}

check_config() {
  local missing=0

  for key in PROJECT_ROOT RAW_FASTQ_DIR SPECIES LIBRARY_TYPES SLURM_ACCOUNT SLURM_QOS; do
    if [ -z "${!key:-}" ]; then
      echo "ERROR: ${key} is not set in ${CONFIG_FILE}" >&2
      missing=1
    fi
  done
  [ "$missing" -eq 0 ] || return 1

  case "$PROJECT_ROOT" in
    /path/to/*) config_die "PROJECT_ROOT is still the placeholder value" ;;
  esac
  case "$RAW_FASTQ_DIR" in
    /path/to/*) config_die "RAW_FASTQ_DIR is still the placeholder value" ;;
  esac
  if [ "$SLURM_ACCOUNT" = "your_slurm_account" ]; then
    config_die "SLURM_ACCOUNT is still the placeholder value"
  fi

  case "$SPECIES" in
    mouse|human) ;;
    *) config_die "SPECIES must be 'mouse' or 'human', got '${SPECIES}'" ;;
  esac

  case " $LIBRARY_TYPES " in
    *" GEX "*) ;;
    *) config_die "LIBRARY_TYPES must include GEX, got '${LIBRARY_TYPES}'" ;;
  esac

  return 0
}

# has_library GEX  ->  returns 0 if that library type is configured
has_library() {
  case " ${LIBRARY_TYPES} " in
    *" $1 "*) return 0 ;;
    *) return 1 ;;
  esac
}

# Strip comments and the header from samples.tsv, emitting only data rows.
samples_data() {
  grep -v '^#' "$SAMPLES_TSV" | awk 'NR>1 && NF>0'
}

# Unique pool names, in the order they first appear in samples.tsv.
sample_pools() {
  samples_data | cut -f1 | awk '!seen[$0]++'
}

# Unique analysis groups, in the order they first appear.
sample_groups() {
  samples_data | cut -f5 | awk '!seen[$0]++'
}

# The fastq prefix for a given pool.
pool_prefix() {
  samples_data | awk -F'\t' -v p="$1" '$1==p {print $2; exit}'
}
