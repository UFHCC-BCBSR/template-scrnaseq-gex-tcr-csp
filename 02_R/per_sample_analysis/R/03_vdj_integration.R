# =============================================================================
# Step 3: Attach V(D)J (TCR / BCR) data to the Seurat object
# =============================================================================
# Summarises the filtered contig annotations per cell barcode and adds the
# result to the object's metadata, so that clonotype and receptor status can be
# used in clustering plots and differential expression.
#
# Skipped gracefully if neither TCR nor BCR is in LIBRARY_TYPES.
#
# Run:
#     Rscript 03_vdj_integration.R [GROUP_NAME]
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tidyr)
})

# Find and load the shared config by walking up to the repo root.
.p <- normalizePath(getwd(), mustWork = FALSE)
while (!file.exists(file.path(.p, "02_R", "R", "config.R")) && dirname(.p) != .p) {
  .p <- dirname(.p)
}
if (!file.exists(file.path(.p, "02_R", "R", "config.R"))) {
  stop("Could not find 02_R/R/config.R. Run this script from inside the repository.")
}
source(file.path(.p, "02_R", "R", "config.R"))

obj      <- readRDS(file.path(OUTPUT_DIR, "objects", "02_filtered.rds"))
vdj_data <- readRDS(file.path(OUTPUT_DIR, "objects", "01_vdj_data.rds"))

message("=== Step 3: V(D)J integration for group '", DATASET, "' ===")

# Always define these, so downstream scripts can rely on the columns existing
# whether or not this experiment has V(D)J libraries.
obj$has_bcr      <- FALSE
obj$has_tcr      <- FALSE
obj$bcr_complete <- FALSE
obj$tcr_complete <- FALSE

if (!has_library("TCR") && !has_library("BCR")) {
  message("No TCR or BCR library configured; skipping V(D)J integration.")
  saveRDS(obj, file.path(OUTPUT_DIR, "objects", "03_with_vdj.rds"))
  message("=== Step 3 complete (nothing to do) ===")
  message("Next: Rscript 04_norm_and_cluster.R ", DATASET)
  quit(save = "no", status = 0)
}


#' Collapse per-contig V(D)J annotations to one row per cell barcode
#' @param vdj_df filtered_contig_annotations.csv contents
#' @param chain_type "BCR" or "TCR"
process_vdj_data <- function(vdj_df, chain_type) {
  if (is.null(vdj_df) || nrow(vdj_df) == 0) {
    message("No ", chain_type, " data found")
    return(NULL)
  }

  message("Processing ", chain_type, ": ", nrow(vdj_df), " contigs")

  # Keep only contigs Cell Ranger is confident about, in real cells, that
  # encode a functional receptor chain.
  vdj_filtered <- vdj_df %>%
    filter(high_confidence == TRUE, productive == TRUE, is_cell == TRUE)

  message("  after filtering: ", nrow(vdj_filtered), " contigs")
  if (nrow(vdj_filtered) == 0) return(NULL)

  if (chain_type == "BCR") {
    vdj_filtered %>%
      group_by(barcode) %>%
      summarise(
        chains            = paste(unique(chain), collapse = ";"),
        heavy             = any(chain == "IGH"),
        light_kappa       = any(chain == "IGK"),
        light_lambda      = any(chain == "IGL"),
        v_genes           = paste(unique(v_gene[!is.na(v_gene)]), collapse = ";"),
        j_genes           = paste(unique(j_gene[!is.na(j_gene)]), collapse = ";"),
        cdr3s             = paste(unique(cdr3[!is.na(cdr3)]), collapse = ";"),
        productive_contigs = n(),
        .groups = "drop"
      ) %>%
      mutate(
        complete = heavy & (light_kappa | light_lambda),
        chain_pairing = case_when(
          light_kappa & !light_lambda ~ "IGH_IGK",
          light_lambda & !light_kappa ~ "IGH_IGL",
          light_kappa & light_lambda  ~ "IGH_IGK_IGL",
          TRUE                        ~ "IGH_only"
        )
      )
  } else {
    vdj_filtered %>%
      group_by(barcode) %>%
      summarise(
        chains            = paste(unique(chain), collapse = ";"),
        alpha             = any(chain == "TRA"),
        beta              = any(chain == "TRB"),
        v_genes           = paste(unique(v_gene[!is.na(v_gene)]), collapse = ";"),
        j_genes           = paste(unique(j_gene[!is.na(j_gene)]), collapse = ";"),
        cdr3s             = paste(unique(cdr3[!is.na(cdr3)]), collapse = ";"),
        productive_contigs = n(),
        .groups = "drop"
      ) %>%
      mutate(
        complete = alpha & beta,
        chain_pairing = case_when(
          alpha & beta   ~ "TRA_TRB",
          alpha & !beta  ~ "TRA_only",
          !alpha & beta  ~ "TRB_only",
          TRUE           ~ "unknown"
        )
      )
  }
}


# Cell names in the merged object are "<sample_id>_<barcode>", while V(D)J
# barcodes are bare. Strip everything up to the final underscore to match.
# Computed once here so it is in scope for both receptor types.
obj_barcode_only <- sub(".*_", "", colnames(obj))

#' Copy a V(D)J summary onto the Seurat object metadata with a prefix
attach_vdj <- function(obj, summary_df, prefix) {
  if (is.null(summary_df)) return(obj)

  idx <- match(obj_barcode_only, summary_df$barcode)

  for (col in setdiff(names(summary_df), "barcode")) {
    obj[[paste0(prefix, "_", col)]] <- summary_df[[col]][idx]
  }

  obj[[paste0("has_", prefix)]] <- !is.na(idx)

  complete_col <- paste0(prefix, "_complete")
  obj[[complete_col]] <- !is.na(idx) &
    !is.na(obj[[complete_col]][, 1]) &
    obj[[complete_col]][, 1]

  message(toupper(prefix), "+ cells: ", sum(!is.na(idx)))
  obj
}

if (has_library("BCR")) {
  obj <- attach_vdj(obj, process_vdj_data(vdj_data$bcr, "BCR"), "bcr")
}
if (has_library("TCR")) {
  obj <- attach_vdj(obj, process_vdj_data(vdj_data$tcr, "TCR"), "tcr")
}

# ---- summary ----
vdj_summary_stats <- obj@meta.data %>%
  group_by(condition, sample_id) %>%
  summarise(
    total_cells  = n(),
    bcr_positive = sum(has_bcr, na.rm = TRUE),
    tcr_positive = sum(has_tcr, na.rm = TRUE),
    bcr_complete = sum(bcr_complete, na.rm = TRUE),
    tcr_complete = sum(tcr_complete, na.rm = TRUE),
    bcr_percent  = round(100 * bcr_positive / total_cells, 1),
    tcr_percent  = round(100 * tcr_positive / total_cells, 1),
    .groups = "drop"
  )
write_csv(vdj_summary_stats, file.path(OUTPUT_DIR, "tables", "03_vdj_summary.csv"))

saveRDS(obj, file.path(OUTPUT_DIR, "objects", "03_with_vdj.rds"))

message("=== Step 3 complete ===")
print(vdj_summary_stats)
message("Next: Rscript 04_norm_and_cluster.R ", DATASET)
