# =============================================================================
# Shared configuration for the R pipeline
# =============================================================================
# Sourced at the top of every R script. Reads config/config.env and
# config/samples.tsv so that no script contains a hard-coded path.
#
# Usage in a script:
#     source(file.path(dirname(dirname(getwd())), "R", "config.R"))
# but in practice just use the one-liner that each script already has, which
# finds this file by walking up from the working directory.
#
# Defines:
#   CFG              named list of everything in config.env
#   REPO_ROOT        absolute path to the repo checkout
#   SAMPLES          data.frame from samples.tsv
#   GROUPS           unique analysis groups
#   DATASET          the group this run is analysing
#   OUTPUT_DIR       where this group's outputs are written
#   WEB_DIR          where results are published
#   WEB_URL          public URL matching WEB_DIR
#   CELLRANGER_OUT   Cell Ranger `outs` dir for this group's pool
#   MT_PATTERN, RIBO_PATTERN, HB_PATTERN   species-appropriate regexes
#   CELL_TYPE_MARKERS                      species-appropriate marker list
#   sample_metadata()                      metadata for one group
#   has_library()                          TRUE if a library type is configured
# =============================================================================


# Small null-coalescing helper: use `b` when `a` is NULL or empty.
`%||%` <- function(a, b) if (is.null(a) || !nzchar(a)) b else a


# ---- locate the repo root ---------------------------------------------------
# Walk up from the working directory until we find config/config.env. This lets
# the scripts be run from anywhere (their own dir, the repo root, or via sbatch).
.find_repo_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = FALSE)
  for (i in 1:10) {
    if (file.exists(file.path(path, "config", "config.env"))) return(path)
    parent <- dirname(path)
    if (parent == path) break
    path <- parent
  }
  stop("Could not locate config/config.env by walking up from '", start,
       "'. Run this script from inside the repository.")
}

REPO_ROOT <- .find_repo_root()


# ---- read config.env --------------------------------------------------------
.read_config_env <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]

  cfg <- list()
  for (line in lines) {
    eq <- regexpr("=", line, fixed = TRUE)
    if (eq < 1) next
    key <- trimws(substr(line, 1, eq - 1))
    val <- trimws(substr(line, eq + 1, nchar(line)))
    # strip one layer of surrounding quotes
    val <- sub('^"(.*)"$', "\\1", val)
    val <- sub("^'(.*)'$", "\\1", val)
    cfg[[key]] <- val
  }
  cfg
}

CFG <- .read_config_env(file.path(REPO_ROOT, "config", "config.env"))


# ---- validate ---------------------------------------------------------------
.required <- c("PROJECT_ROOT", "SPECIES", "LIBRARY_TYPES", "PROJECT_NAME")
.missing <- .required[!.required %in% names(CFG)]
if (length(.missing) > 0) {
  stop("config.env is missing required keys: ", paste(.missing, collapse = ", "))
}
if (grepl("^/path/to/", CFG$PROJECT_ROOT)) {
  stop("PROJECT_ROOT in config/config.env is still the placeholder value. ",
       "Edit config/config.env before running the pipeline.")
}
if (!CFG$SPECIES %in% c("mouse", "human")) {
  stop("SPECIES must be 'mouse' or 'human', got '", CFG$SPECIES, "'")
}


# ---- library types ----------------------------------------------------------
LIBRARY_TYPES <- strsplit(trimws(CFG$LIBRARY_TYPES), "\\s+")[[1]]

#' Is a given library type part of this experiment?
#' @param lib one of "GEX", "CSP", "TCR", "BCR"
has_library <- function(lib) lib %in% LIBRARY_TYPES


# ---- sample sheet -----------------------------------------------------------
.samples_path <- file.path(REPO_ROOT, "config", "samples.tsv")
SAMPLES <- read.delim(
  .samples_path,
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE,
  colClasses = "character"
)

.expected_cols <- c("pool", "pool_prefix", "hashtag_id", "sample_id",
                    "group", "condition", "replicate")
.missing_cols <- setdiff(.expected_cols, colnames(SAMPLES))
if (length(.missing_cols) > 0) {
  stop("config/samples.tsv is missing columns: ", paste(.missing_cols, collapse = ", "))
}
if (nrow(SAMPLES) == 0) {
  stop("config/samples.tsv contains no data rows.")
}
if (any(duplicated(SAMPLES$sample_id))) {
  stop("config/samples.tsv has duplicate sample_id values: ",
       paste(unique(SAMPLES$sample_id[duplicated(SAMPLES$sample_id)]), collapse = ", "))
}

GROUPS <- unique(SAMPLES$group)
POOLS  <- unique(SAMPLES$pool)


# ---- which group are we analysing? ------------------------------------------
# Priority: explicit DATASET already set by the caller > command-line argument >
# DATASET environment variable > first group in the sample sheet.
if (!exists("DATASET") || is.null(DATASET) || !nzchar(DATASET)) {
  .args <- commandArgs(trailingOnly = TRUE)
  .env_dataset <- Sys.getenv("DATASET", unset = "")
  DATASET <- if (length(.args) > 0 && nzchar(.args[1])) {
    .args[1]
  } else if (nzchar(.env_dataset)) {
    .env_dataset
  } else {
    GROUPS[1]
  }
}

if (!DATASET %in% GROUPS) {
  stop("Group '", DATASET, "' is not present in config/samples.tsv.\n",
       "  Available groups: ", paste(GROUPS, collapse = ", "), "\n",
       "  Pass one as an argument, e.g.  Rscript 01_load_data.R ", GROUPS[1])
}


# ---- derived paths ----------------------------------------------------------
PROJECT_ROOT <- CFG$PROJECT_ROOT

# Cell Ranger writes to <pool>_run/<pool>_multi/outs inside the checkout.
# A group maps to the pool(s) its samples came from.
.group_pools <- unique(SAMPLES$pool[SAMPLES$group == DATASET])
CELLRANGER_BASE <- file.path(REPO_ROOT, "01_preprocessing", "02_cellranger")

#' Cell Ranger `outs` directory for a pool
cellranger_outs <- function(pool) {
  file.path(CELLRANGER_BASE, paste0(pool, "_run"), paste0(pool, "_multi"), "outs")
}

# The primary pool for this group (most groups map to exactly one).
CELLRANGER_OUT <- cellranger_outs(.group_pools[1])

ANALYSIS_DIR <- file.path(REPO_ROOT, "02_R", "per_sample_analysis")
OUTPUT_DIR   <- file.path(ANALYSIS_DIR, "analysis_outputs", DATASET)

WEB_DIR <- file.path(CFG$WEB_PUBLISH_DIR %||% "", DATASET)
WEB_URL <- paste0(sub("/$", "", CFG$WEB_PUBLISH_URL %||% ""), "/", DATASET)


# ---- species-specific settings ----------------------------------------------
# Mouse gene symbols are Xkr4-style (first letter capitalised); human are
# XKR4-style (all caps). Both the QC regexes and the marker lists differ.
if (CFG$SPECIES == "mouse") {
  MT_PATTERN   <- "^mt-"
  RIBO_PATTERN <- "^Rpl|^Rps"
  HB_PATTERN   <- "^Hb[^(p)]"

  CELL_TYPE_MARKERS <- list(
    "T_cells"         = c("Cd3e", "Cd3d", "Cd3g"),
    "CD4_T"           = c("Cd4"),
    "CD8_T"           = c("Cd8a", "Cd8b1"),
    "Regulatory_T"    = c("Foxp3", "Il2ra", "Ctla4", "Ikzf2"),
    "NK_cells"        = c("Nkg7", "Gzma", "Klrb1c", "Ncr1"),
    "B_cells"         = c("Cd19", "Ms4a1", "Cd79a", "Pax5"),
    "Plasma_cells"    = c("Jchain", "Mzb1", "Xbp1", "Sdc1"),
    "Macrophages"     = c("Adgre1", "Cd68", "Itgam", "Fcgr1"),
    "Dendritic_cells" = c("Itgax", "Flt3", "Zbtb46"),
    "Monocytes"       = c("Ly6c2", "Ccr2", "Nr4a1"),
    "Neutrophils"     = c("S100a8", "S100a9", "Ly6g")
  )
} else {
  MT_PATTERN   <- "^MT-"
  RIBO_PATTERN <- "^RPL|^RPS"
  HB_PATTERN   <- "^HB[^(P)]"

  CELL_TYPE_MARKERS <- list(
    "T_cells"         = c("CD3E", "CD3D", "CD3G"),
    "CD4_T"           = c("CD4", "IL7R"),
    "CD8_T"           = c("CD8A", "CD8B"),
    "Regulatory_T"    = c("FOXP3", "IL2RA", "CTLA4", "IKZF2"),
    "NK_cells"        = c("NKG7", "GNLY", "KLRD1", "NCAM1"),
    "B_cells"         = c("CD19", "MS4A1", "CD79A", "PAX5"),
    "Plasma_cells"    = c("JCHAIN", "MZB1", "XBP1", "SDC1"),
    "Macrophages"     = c("CD68", "ITGAM", "FCGR1A", "MRC1"),
    "Dendritic_cells" = c("ITGAX", "FLT3", "ZBTB46"),
    "Monocytes"       = c("CD14", "FCGR3A", "CCR2"),
    "Neutrophils"     = c("S100A8", "S100A9", "FCGR3B")
  )
}


# ---- helpers ----------------------------------------------------------------

#' Metadata for the samples in one analysis group
#' @param group group name; defaults to the current DATASET
#' @return data.frame with one row per sample
sample_metadata <- function(group = DATASET) {
  md <- SAMPLES[SAMPLES$group == group, , drop = FALSE]
  if (nrow(md) == 0) {
    stop("No samples found for group '", group, "' in config/samples.tsv")
  }
  data.frame(
    sample_id = md$sample_id,
    condition = md$condition,
    replicate = md$replicate,
    group     = md$group,
    pool      = md$pool,
    hashtag   = md$hashtag_id,
    stringsAsFactors = FALSE
  )
}

#' Create the standard output subdirectories for a group
init_output_dirs <- function(dir = OUTPUT_DIR) {
  for (sub in c("plots", "tables", "objects")) {
    dir.create(file.path(dir, sub), recursive = TRUE, showWarnings = FALSE)
  }
  invisible(dir)
}


message("--- config loaded ---")
message("  repo:     ", REPO_ROOT)
message("  group:    ", DATASET, "  (", nrow(sample_metadata()), " samples)")
message("  species:  ", CFG$SPECIES)
message("  libs:     ", paste(LIBRARY_TYPES, collapse = ", "))
message("  outputs:  ", OUTPUT_DIR)
