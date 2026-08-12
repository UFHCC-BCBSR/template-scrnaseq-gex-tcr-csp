# =============================================================================
# Integration: export the integrated dataset to Loupe Browser format
# =============================================================================
# Writes a single .cloupe file spanning all analysis groups, so collaborators
# can explore the combined dataset interactively without R.
#
# Requires the loupeR package: https://github.com/10XGenomics/loupeR
#
# Run:
#     Rscript 02_loupe.R
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

OUT <- file.path(REPO_ROOT, "02_R", "integrated_analysis", "analysis_outputs")
obj_path <- file.path(OUT, "objects", "integrated.rds")

if (!file.exists(obj_path)) {
  stop("Integrated object not found:\n  ", obj_path,
       "\nRun 01_integrate_data.R first.")
}

message("=== Loupe export: integrated dataset ===")
obj <- readRDS(obj_path)

message("Assays: ",     paste(names(obj@assays), collapse = ", "))
message("Reductions: ", paste(names(obj@reductions), collapse = ", "))

DefaultAssay(obj) <- "RNA"

# Seurat v5 keeps per-sample layers after a merge; Loupe needs one matrix.
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

# Loupe requires character metadata, not factors.
for (col in colnames(obj@meta.data)) {
  if (is.factor(obj@meta.data[[col]])) {
    obj@meta.data[[col]] <- as.character(obj@meta.data[[col]])
  }
}

# Categorical columns to filter on inside Loupe.
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

# Combined group and condition label, useful for side-by-side comparison.
if (all(c("analysis_group", "condition") %in% colnames(obj@meta.data))) {
  obj$group_condition <- paste(obj$analysis_group, obj$condition, sep = "_")
}

output_file <- file.path(OUT, "objects", paste0(CFG$PROJECT_NAME, "_integrated.cloupe"))
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

message("Creating Loupe file...")
create_loupe_from_seurat(obj, output_name = output_file, force = TRUE)

if (!file.exists(output_file)) output_file <- paste0(output_file, ".cloupe")

message("Wrote: ", output_file)
message("Size: ", round(file.size(output_file) / 1024^2, 2), " MB")

# Publish alongside the per-group results, if a web directory is configured.
if (!grepl("^/path/to/", CFG$WEB_PUBLISH_DIR %||% "/path/to/")) {
  web_integrated <- file.path(CFG$WEB_PUBLISH_DIR, "integrated")
  dir.create(web_integrated, recursive = TRUE, showWarnings = FALSE)
  file.copy(output_file, file.path(web_integrated, basename(output_file)), overwrite = TRUE)
  message("Published to: ", file.path(web_integrated, basename(output_file)))
}

group_summary <- obj@meta.data %>%
  group_by(analysis_group) %>%
  summarise(
    cells   = n(),
    tcr_pos = sum(has_tcr %in% TRUE),
    bcr_pos = sum(has_bcr %in% TRUE),
    .groups = "drop"
  )
write.csv(group_summary,
          file.path(OUT, "tables", "loupe_group_breakdown.csv"),
          row.names = FALSE)

message("\n=== Loupe export complete ===")
print(group_summary)
