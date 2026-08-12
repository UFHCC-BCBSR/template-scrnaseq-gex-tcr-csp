# =============================================================================
# Integration: combine all analysis groups with Harmony batch correction
# =============================================================================
# Merges the per-group Seurat objects and corrects for group-level batch
# effects using Harmony, so that clusters reflect cell identity rather than
# which group a cell came from.
#
# Requires every group to have completed at least step 3 of the per-sample
# pipeline.
#
# Run:
#     Rscript 01_integrate_data.R
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(readr)
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


# ---- parameters -------------------------------------------------------------
N_VARIABLE_FEATURES <- 2000
N_DIMS              <- 30
CLUSTER_RESOLUTION  <- 0.5


PER_SAMPLE_OUTPUTS <- file.path(REPO_ROOT, "02_R", "per_sample_analysis", "analysis_outputs")
OUT <- file.path(REPO_ROOT, "02_R", "integrated_analysis", "analysis_outputs")
for (sub in c("objects", "plots", "tables")) {
  dir.create(file.path(OUT, sub), recursive = TRUE, showWarnings = FALSE)
}

message("=== Integrating ", length(GROUPS), " groups ===")
message("Groups: ", paste(GROUPS, collapse = ", "))

if (length(GROUPS) < 2) {
  stop("Only one group is defined in config/samples.tsv, so there is nothing ",
       "to integrate. This step is for combining two or more groups.")
}

# ---- load each group's object ----
seurat_objects <- list()
vdj_summaries  <- list()

for (group in GROUPS) {
  obj_path <- file.path(PER_SAMPLE_OUTPUTS, group, "objects", "03_with_vdj.rds")
  if (!file.exists(obj_path)) {
    stop("Missing input for group '", group, "':\n  ", obj_path,
         "\nRun the per-sample pipeline for that group first.")
  }

  message("Loading ", group, "...")
  o <- readRDS(obj_path)
  o$analysis_group <- group
  seurat_objects[[group]] <- o

  vdj_path <- file.path(PER_SAMPLE_OUTPUTS, group, "tables", "03_vdj_summary.csv")
  if (file.exists(vdj_path)) {
    vdj_summaries[[group]] <- read_csv(vdj_path, show_col_types = FALSE)
  }

  message("  cells: ", ncol(o),
          "  TCR+: ", sum(o$has_tcr, na.rm = TRUE),
          "  BCR+: ", sum(o$has_bcr, na.rm = TRUE))
}

# ---- merge and preprocess ----
message("\nMerging...")
integrated_obj <- merge(
  seurat_objects[[1]],
  y = seurat_objects[-1],
  add.cell.ids = names(seurat_objects)
)

message("Normalising and running PCA...")
integrated_obj <- NormalizeData(integrated_obj, normalization.method = "LogNormalize",
                                scale.factor = 10000)
integrated_obj <- FindVariableFeatures(integrated_obj, selection.method = "vst",
                                       nfeatures = N_VARIABLE_FEATURES)
integrated_obj <- ScaleData(integrated_obj, verbose = FALSE)
integrated_obj <- RunPCA(integrated_obj, npcs = 50, verbose = FALSE)

pdf(file.path(OUT, "plots", "integration_elbow_plot.pdf"))
print(ElbowPlot(integrated_obj, ndims = 50))
dev.off()

# ---- Harmony ----
# Corrects the PCA embedding so that cells cluster by biology rather than by
# which group they came from. Compare the before and after UMAPs to check it
# worked: groups should overlap afterwards, not sit in separate islands.
message("Running Harmony on 'analysis_group'...")
integrated_obj <- RunHarmony(
  integrated_obj,
  group.by.vars = "analysis_group",
  assay.use     = "RNA",
  verbose       = TRUE
)

message("Clustering and running UMAP...")
integrated_obj <- RunUMAP(integrated_obj, reduction = "harmony", dims = 1:N_DIMS)
integrated_obj <- FindNeighbors(integrated_obj, reduction = "harmony", dims = 1:N_DIMS)
integrated_obj <- FindClusters(integrated_obj, resolution = CLUSTER_RESOLUTION)

integrated_obj <- JoinLayers(integrated_obj)

# ---- plots ----
message("Writing plots...")

p1 <- DimPlot(integrated_obj, reduction = "umap", group.by = "analysis_group") +
  ggtitle("After Harmony: cells by group")
ggsave(file.path(OUT, "plots", "integration_by_group.pdf"), p1, width = 10, height = 8)

p2 <- DimPlot(integrated_obj, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, repel = TRUE) + ggtitle("Integrated clusters")
ggsave(file.path(OUT, "plots", "integration_clusters.pdf"), p2, width = 10, height = 8)

if ("condition" %in% colnames(integrated_obj@meta.data)) {
  p3 <- DimPlot(integrated_obj, reduction = "umap", group.by = "condition") +
    ggtitle("Cells by condition")
  ggsave(file.path(OUT, "plots", "integration_by_condition.pdf"), p3, width = 10, height = 8)
}

vdj_plots <- list()
if (sum(integrated_obj$has_tcr, na.rm = TRUE) > 0) {
  vdj_plots$tcr <- DimPlot(integrated_obj, reduction = "umap", group.by = "has_tcr",
                           cols = c("grey90", "steelblue")) + ggtitle("TCR+ cells")
}
if (sum(integrated_obj$has_bcr, na.rm = TRUE) > 0) {
  vdj_plots$bcr <- DimPlot(integrated_obj, reduction = "umap", group.by = "has_bcr",
                           cols = c("grey90", "firebrick")) + ggtitle("BCR+ cells")
}
if (length(vdj_plots) > 0) {
  ggsave(file.path(OUT, "plots", "integration_vdj_overview.pdf"),
         patchwork::wrap_plots(vdj_plots, ncol = length(vdj_plots)),
         width = 8 * length(vdj_plots), height = 6)
}

# ---- tables ----
integration_summary <- integrated_obj@meta.data %>%
  group_by(analysis_group) %>%
  summarise(
    total_cells  = n(),
    tcr_positive = sum(has_tcr, na.rm = TRUE),
    bcr_positive = sum(has_bcr, na.rm = TRUE),
    tcr_percent  = round(100 * tcr_positive / total_cells, 1),
    bcr_percent  = round(100 * bcr_positive / total_cells, 1),
    .groups = "drop"
  )
write_csv(integration_summary, file.path(OUT, "tables", "integration_summary.csv"))

# Cluster by group: a cluster made up almost entirely of one group may be a
# residual batch effect rather than a real population.
cluster_by_group <- as.data.frame.matrix(
  table(integrated_obj$seurat_clusters, integrated_obj$analysis_group)
)
cluster_by_group$cluster <- rownames(cluster_by_group)
write_csv(cluster_by_group, file.path(OUT, "tables", "cluster_by_group.csv"))

if (length(vdj_summaries) > 0) {
  write_csv(bind_rows(vdj_summaries, .id = "analysis_group"),
            file.path(OUT, "tables", "combined_vdj_summary.csv"))
}

saveRDS(integrated_obj, file.path(OUT, "objects", "integrated.rds"))

message("\n=== Integration complete ===")
message("Total cells: ", ncol(integrated_obj))
message("Clusters: ", length(levels(integrated_obj$seurat_clusters)))
print(integration_summary)
message("\nCheck integration_by_group.pdf: groups should overlap, not form ",
        "separate islands.")
message("Next: Rscript 02_loupe.R")
