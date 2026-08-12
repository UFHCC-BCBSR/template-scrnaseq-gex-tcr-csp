# =============================================================================
# Step 9: Export to Loupe Browser format
# =============================================================================
# Writes a .cloupe file so collaborators without R can explore the data
# interactively in 10x Loupe Browser: filter cells by metadata, view clusters,
# and look up gene expression.
#
# Requires the loupeR package: https://github.com/10XGenomics/loupeR
#
# Run:
#     Rscript 09_loupe.R [GROUP_NAME]
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
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

if (!requireNamespace("loupeR", quietly = TRUE)) {
  stop("The loupeR package is not installed.\n",
       "Install it from https://github.com/10XGenomics/loupeR, or skip this step.")
}
library(loupeR)

message("=== Step 9: Loupe export for group '", DATASET, "' ===")

obj <- readRDS(file.path(OUTPUT_DIR, "objects", "06_final_analysis.rds"))

DefaultAssay(obj) <- "RNA"

# Seurat v5 keeps per-sample layers after a merge; Loupe needs one matrix.
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

# Loupe requires character metadata, not factors.
for (col in colnames(obj@meta.data)) {
  if (is.factor(obj@meta.data[[col]])) {
    obj@meta.data[[col]] <- as.character(obj@meta.data[[col]])
  }
}

# Add readable categorical columns to filter on inside Loupe.
obj$cluster_label <- paste0("Cluster_", obj$seurat_clusters)

if (has_library("TCR")) {
  obj$tcr_status <- ifelse(obj$has_tcr %in% TRUE, "TCR_positive", "TCR_negative")
}
if (has_library("BCR")) {
  obj$bcr_status <- ifelse(obj$has_bcr %in% TRUE, "BCR_positive", "BCR_negative")
}
if (has_library("TCR") || has_library("BCR")) {
  obj$receptor_type <- dplyr::case_when(
    obj$has_bcr %in% TRUE ~ "B_cell_receptor",
    obj$has_tcr %in% TRUE ~ "T_cell_receptor",
    TRUE                  ~ "None_detected"
  )
}

# Write next to the other outputs, and to the web directory if configured.
out_dir <- file.path(OUTPUT_DIR, "objects")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(out_dir, paste0(CFG$PROJECT_NAME, "_", DATASET, ".cloupe"))

message("Creating Loupe file...")
create_loupe_from_seurat(obj, output_name = output_file, force = TRUE)

# create_loupe_from_seurat appends .cloupe if it is not already there.
if (!file.exists(output_file)) output_file <- paste0(output_file, ".cloupe")

message("Wrote: ", output_file)
message("Size: ", round(file.size(output_file) / 1024^2, 2), " MB")

if (!grepl("^/path/to/", CFG$WEB_PUBLISH_DIR %||% "/path/to/")) {
  dir.create(WEB_DIR, recursive = TRUE, showWarnings = FALSE)
  file.copy(output_file, file.path(WEB_DIR, basename(output_file)), overwrite = TRUE)
  message("Published to: ", file.path(WEB_DIR, basename(output_file)))
}

summary_info <- data.frame(
  group        = DATASET,
  total_cells  = ncol(obj),
  total_genes  = nrow(obj),
  clusters     = length(unique(obj$seurat_clusters)),
  tcr_positive = sum(obj$has_tcr %in% TRUE),
  bcr_positive = sum(obj$has_bcr %in% TRUE),
  file_size_mb = round(file.size(output_file) / 1024^2, 2)
)
write.csv(summary_info,
          file.path(OUTPUT_DIR, "tables", "09_loupe_summary.csv"),
          row.names = FALSE)

message("=== Step 9 complete ===")
print(summary_info)
