# =============================================================================
# Step 1: Load Cell Ranger output into a Seurat object
# =============================================================================
# Reads the per-sample matrices produced by cellranger multi for one analysis
# group, attaches hashtag and V(D)J data, merges the samples, and saves the
# result for the next step.
#
# Run:
#     Rscript 01_load_data.R              # first group in config/samples.tsv
#     Rscript 01_load_data.R GROUP_NAME   # a specific group
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(stringr)
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

init_output_dirs()

metadata <- sample_metadata()
samples  <- metadata$sample_id

message("=== Step 1: Loading data for group '", DATASET, "' ===")
message("Cell Ranger output: ", CELLRANGER_OUT)

if (!dir.exists(CELLRANGER_OUT)) {
  stop("Cell Ranger output not found:\n  ", CELLRANGER_OUT,
       "\nHas the preprocessing finished? See 01_preprocessing/02_cellranger/")
}


#' Load one demultiplexed sample: counts, hashtags and V(D)J contigs
load_sample_data <- function(sample_id, base_path) {
  message("Loading sample: ", sample_id)

  sample_dir  <- file.path(base_path, "per_sample_outs", sample_id)
  count_path  <- file.path(sample_dir, "count", "sample_filtered_feature_bc_matrix")
  vdj_b_path  <- file.path(sample_dir, "vdj_b", "filtered_contig_annotations.csv")
  vdj_t_path  <- file.path(sample_dir, "vdj_t", "filtered_contig_annotations.csv")
  metrics_path <- file.path(sample_dir, "metrics_summary.csv")

  if (!dir.exists(count_path)) {
    stop("Count matrix not found for sample '", sample_id, "':\n  ", count_path,
         "\nCheck that sample_id in config/samples.tsv matches the Cell Ranger ",
         "per_sample_outs directory names.")
  }

  counts <- Read10X(count_path)

  # Read10X returns a list when the library has multiple feature types
  # (Gene Expression + Antibody Capture), or a plain matrix when it is GEX only.
  if (is.list(counts)) {
    gex_counts <- counts[["Gene Expression"]]
    hto_counts <- if ("Antibody Capture" %in% names(counts)) counts[["Antibody Capture"]] else NULL
  } else {
    gex_counts <- counts
    hto_counts <- NULL
  }

  # Keep only barcodes seen in both assays, so the matrices stay aligned.
  if (!is.null(hto_counts)) {
    joint_bcs <- intersect(colnames(gex_counts), colnames(hto_counts))
    message("  - Gene Expression: ", ncol(gex_counts), " barcodes")
    message("  - Antibody Capture: ", ncol(hto_counts), " barcodes")
    message("  - Joint barcodes: ", length(joint_bcs))
    gex_counts <- gex_counts[, joint_bcs, drop = FALSE]
    hto_counts <- hto_counts[, joint_bcs, drop = FALSE]
  }

  obj <- CreateSeuratObject(
    counts     = gex_counts,
    project    = sample_id,
    min.cells  = 3,
    min.features = 200
  )

  # Attach hashtag counts, subset to the cells that survived the filter above.
  if (!is.null(hto_counts)) {
    hto_filtered <- hto_counts[, colnames(obj), drop = FALSE]
    obj[["HTO"]] <- CreateAssayObject(counts = hto_filtered)
    message("  - Added HTO assay: ", nrow(hto_filtered), " features, ",
            ncol(hto_filtered), " cells")
  }

  # Attach this sample's metadata to every one of its cells.
  info <- metadata[metadata$sample_id == sample_id, , drop = FALSE]
  for (col in setdiff(names(info), "sample_id")) {
    obj[[col]] <- info[[col]]
  }
  obj[["sample_id"]] <- sample_id

  # V(D)J contigs, if those libraries were sequenced.
  vdj_data <- list()
  if (has_library("BCR") && file.exists(vdj_b_path)) {
    vdj_data$bcr <- read_csv(vdj_b_path, show_col_types = FALSE)
    message("  - BCR contigs: ", nrow(vdj_data$bcr))
  }
  if (has_library("TCR") && file.exists(vdj_t_path)) {
    vdj_data$tcr <- read_csv(vdj_t_path, show_col_types = FALSE)
    message("  - TCR contigs: ", nrow(vdj_data$tcr))
  }
  if (file.exists(metrics_path)) {
    vdj_data$metrics <- read_csv(metrics_path, show_col_types = FALSE)
  }

  list(seurat = obj, vdj = vdj_data)
}


# ---- load every sample in the group ----
sample_data <- list()
for (sid in samples) {
  sample_data[[sid]] <- load_sample_data(sid, CELLRANGER_OUT)
}

seurat_objects <- lapply(sample_data, function(x) x$seurat)

message("Merging ", length(seurat_objects), " samples...")
merged_obj <- merge(
  seurat_objects[[1]],
  y = if (length(seurat_objects) > 1) seurat_objects[-1] else NULL,
  add.cell.ids = names(seurat_objects),
  project = paste0(CFG$PROJECT_NAME, "_", DATASET)
)

message("Merged object: ", ncol(merged_obj), " cells, ", nrow(merged_obj), " genes")

# NOTE: sample_id is set per-object before merging (above) rather than parsed
# back out of the merged cell names. Cell names look like <sample_id>_<barcode>,
# but sample ids can themselves contain underscores, so parsing them out with a
# regex silently breaks for many naming schemes.

# ---- combine V(D)J across samples ----
message("Compiling V(D)J data...")
all_vdj_data <- list(
  bcr     = if (has_library("BCR")) bind_rows(lapply(sample_data, function(x) x$vdj$bcr)) else NULL,
  tcr     = if (has_library("TCR")) bind_rows(lapply(sample_data, function(x) x$vdj$tcr)) else NULL,
  metrics = bind_rows(lapply(sample_data, function(x) x$vdj$metrics))
)

# ---- save ----
saveRDS(merged_obj,  file.path(OUTPUT_DIR, "objects", "01_merged_raw.rds"))
saveRDS(all_vdj_data, file.path(OUTPUT_DIR, "objects", "01_vdj_data.rds"))
saveRDS(metadata,    file.path(OUTPUT_DIR, "objects", "01_sample_metadata.rds"))

initial_summary <- data.frame(
  sample_id = samples,
  n_cells   = sapply(seurat_objects, ncol),
  n_genes   = sapply(seurat_objects, nrow),
  condition = metadata$condition,
  replicate = metadata$replicate
)
write_csv(initial_summary, file.path(OUTPUT_DIR, "tables", "01_sample_summary.csv"))

message("=== Step 1 complete ===")
print(initial_summary)
message("Next: Rscript 02_qc_filtering.R ", DATASET)
