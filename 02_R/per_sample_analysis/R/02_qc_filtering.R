# =============================================================================
# Step 2: Quality control and cell filtering
# =============================================================================
# Computes per-cell QC metrics, plots them, and removes low-quality cells.
#
# READ THE PLOTS BEFORE TRUSTING THE THRESHOLDS. The defaults below are
# reasonable starting points, not universal truths. Different tissues and
# protocols need different cutoffs.
#
# Run:
#     Rscript 02_qc_filtering.R [GROUP_NAME]
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(patchwork)
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


# ---- QC thresholds ----------------------------------------------------------
# EDIT THESE after looking at the "before filtering" plots this script writes.
MIN_GENES <- 200    # cells with fewer detected genes are usually debris
MAX_GENES <- 5000   # cells with far more may be doublets
MIN_UMIS  <- 1000   # minimum total counts
MAX_MT    <- 20     # max % mitochondrial reads; dying cells run high
MAX_HB    <- 5      # max % haemoglobin; flags red blood cell contamination


obj <- readRDS(file.path(OUTPUT_DIR, "objects", "01_merged_raw.rds"))
metadata <- readRDS(file.path(OUTPUT_DIR, "objects", "01_sample_metadata.rds"))

message("=== Step 2: QC and filtering for group '", DATASET, "' ===")
message("Starting with ", ncol(obj), " cells")

# ---- compute QC metrics ----
# Patterns come from config.R and follow the species' gene naming convention.
obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = MT_PATTERN)
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = RIBO_PATTERN)
obj[["percent.hb"]]   <- PercentageFeatureSet(obj, pattern = HB_PATTERN)

if (all(obj[["percent.mt"]] == 0)) {
  warning("No mitochondrial genes matched the pattern '", MT_PATTERN, "'. ",
          "Check that SPECIES in config/config.env matches your reference.")
}

# ---- summarise before filtering ----
qc_summary <- obj@meta.data %>%
  group_by(sample_id, condition) %>%
  summarise(
    n_cells      = n(),
    median_genes = median(nFeature_RNA),
    median_umis  = median(nCount_RNA),
    median_mt    = median(percent.mt),
    median_ribo  = median(percent.ribo),
    .groups = "drop"
  )
write_csv(qc_summary, file.path(OUTPUT_DIR, "tables", "02_qc_summary_before.csv"))
print(qc_summary)

# ---- plots ----
p1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = "sample_id", ncol = 3, pt.size = 0) +
  plot_annotation(title = "QC metrics by sample (before filtering)")

p2 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = "condition", ncol = 3, pt.size = 0.1) +
  plot_annotation(title = "QC metrics by condition (before filtering)")

p3 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt",
                     group.by = "condition") + ggtitle("UMIs vs mitochondrial %")
p4 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                     group.by = "condition") + ggtitle("UMIs vs genes detected")

ggsave(file.path(OUTPUT_DIR, "plots", "02_qc_by_sample.pdf"),    p1, width = 15, height = 10)
ggsave(file.path(OUTPUT_DIR, "plots", "02_qc_by_condition.pdf"), p2, width = 12, height = 8)
ggsave(file.path(OUTPUT_DIR, "plots", "02_qc_scatter.pdf"), (p3 | p4), width = 15, height = 6)

# ---- filter ----
message("\nApplying filters:")
message("  ", MIN_GENES, " < nFeature_RNA < ", MAX_GENES)
message("  nCount_RNA > ", MIN_UMIS)
message("  percent.mt < ", MAX_MT)
message("  percent.hb < ", MAX_HB)

obj_filtered <- subset(
  obj,
  subset = nFeature_RNA > MIN_GENES &
           nFeature_RNA < MAX_GENES &
           nCount_RNA   > MIN_UMIS &
           percent.mt   < MAX_MT &
           percent.hb   < MAX_HB
)

message("After filtering: ", ncol(obj_filtered), " cells retained (",
        round(100 * ncol(obj_filtered) / ncol(obj), 1), "%)")

# ---- per-sample retention ----
filter_summary <- data.frame(
  sample_id    = metadata$sample_id,
  condition    = metadata$condition,
  cells_before = sapply(metadata$sample_id, function(x) sum(obj$sample_id == x, na.rm = TRUE)),
  cells_after  = sapply(metadata$sample_id, function(x) sum(obj_filtered$sample_id == x, na.rm = TRUE))
) %>%
  mutate(
    cells_removed    = cells_before - cells_after,
    percent_retained = round(100 * cells_after / cells_before, 2)
  )
write_csv(filter_summary, file.path(OUTPUT_DIR, "tables", "02_filter_summary.csv"))

# A sample losing far more cells than its peers usually means a technical
# problem with that specimen, not biology. Worth investigating before moving on.
if (any(filter_summary$percent_retained < 50, na.rm = TRUE)) {
  warning("Some samples retained under 50% of cells. Review 02_filter_summary.csv ",
          "and the QC plots before continuing.")
}

p5 <- VlnPlot(obj_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = "condition", ncol = 3, pt.size = 0.1) +
  plot_annotation(title = "QC metrics after filtering")
ggsave(file.path(OUTPUT_DIR, "plots", "02_qc_after_filtering.pdf"), p5, width = 12, height = 8)

saveRDS(obj_filtered, file.path(OUTPUT_DIR, "objects", "02_filtered.rds"))

message("=== Step 2 complete ===")
print(filter_summary)
message("Next: Rscript 03_vdj_integration.R ", DATASET)
