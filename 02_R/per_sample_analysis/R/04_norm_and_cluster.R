# =============================================================================
# Step 4: Normalisation, dimensionality reduction and clustering
# =============================================================================
# Standard Seurat workflow: normalise, find variable genes, scale, PCA, then
# neighbour graph, clustering and UMAP.
#
# CHECK THE ELBOW PLOT. N_DIMS below should cover the PCs that carry real
# signal. The default of 25 is a common choice but is not derived from your
# data. Re-run this step after looking at 04_elbow_plot.pdf.
#
# Run:
#     Rscript 04_norm_and_cluster.R [GROUP_NAME]
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
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


# ---- parameters to tune -----------------------------------------------------
N_VARIABLE_FEATURES <- 2000  # genes carried into PCA
N_DIMS              <- 25    # PCs used for clustering and UMAP; see elbow plot
CLUSTER_RESOLUTION  <- 0.5   # higher = more, smaller clusters

set.seed(42)  # UMAP and clustering are stochastic; fix the seed so reruns match


obj <- readRDS(file.path(OUTPUT_DIR, "objects", "03_with_vdj.rds"))

message("=== Step 4: Normalisation and clustering for group '", DATASET, "' ===")

message("Normalising...")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = N_VARIABLE_FEATURES)

top10 <- head(VariableFeatures(obj), 10)
vf_plot <- LabelPoints(plot = VariableFeaturePlot(obj), points = top10, repel = TRUE)
ggsave(file.path(OUTPUT_DIR, "plots", "04_variable_features.pdf"), vf_plot, width = 12, height = 8)

message("Scaling and running PCA...")
obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj))

# The elbow plot tells you how many PCs are worth keeping: look for where the
# curve flattens. If it flattens well before or after N_DIMS, change N_DIMS.
elbow_plot <- ElbowPlot(obj, ndims = 50)
ggsave(file.path(OUTPUT_DIR, "plots", "04_elbow_plot.pdf"), elbow_plot, width = 10, height = 6)

pca_plots <- (DimPlot(obj, reduction = "pca", group.by = "condition") + ggtitle("Condition") |
              DimPlot(obj, reduction = "pca", group.by = "sample_id") + ggtitle("Sample"))
ggsave(file.path(OUTPUT_DIR, "plots", "04_pca_overview.pdf"), pca_plots, width = 15, height = 6)

message("Clustering and running UMAP with ", N_DIMS, " PCs...")
obj <- FindNeighbors(obj, dims = 1:N_DIMS)
obj <- FindClusters(obj, resolution = CLUSTER_RESOLUTION)
obj <- RunUMAP(obj, dims = 1:N_DIMS)

# Seurat v5 splits counts into per-sample layers on merge; rejoin them so that
# downstream differential expression sees a single matrix.
obj <- JoinLayers(obj)

# ---- overview plots ----
p1 <- DimPlot(obj, reduction = "umap", label = TRUE, label.size = 6) +
  ggtitle("Clusters") + NoLegend()
p2 <- DimPlot(obj, reduction = "umap", group.by = "condition") + ggtitle("Condition")
p3 <- DimPlot(obj, reduction = "umap", group.by = "sample_id") + ggtitle("Sample")

# If one sample occupies its own cluster, suspect a batch effect rather than
# biology. The integrated analysis applies Harmony to correct for this.
umap_overview <- (p1 + p2) / (p3 + plot_spacer())
ggsave(file.path(OUTPUT_DIR, "plots", "04_umap_overview.pdf"), umap_overview, width = 15, height = 12)

vdj_plots <- list()
if (sum(obj$has_tcr, na.rm = TRUE) > 0) {
  vdj_plots$tcr <- DimPlot(obj, reduction = "umap", group.by = "has_tcr") + ggtitle("TCR+ cells")
}
if (sum(obj$has_bcr, na.rm = TRUE) > 0) {
  vdj_plots$bcr <- DimPlot(obj, reduction = "umap", group.by = "has_bcr") + ggtitle("BCR+ cells")
}
if (length(vdj_plots) > 0) {
  ggsave(file.path(OUTPUT_DIR, "plots", "04_umap_vdj.pdf"),
         wrap_plots(vdj_plots, ncol = length(vdj_plots)),
         width = 8 * length(vdj_plots), height = 6)
}

split_plot <- DimPlot(obj, reduction = "umap", split.by = "condition") +
  ggtitle("Cells by condition")
ggsave(file.path(OUTPUT_DIR, "plots", "04_umap_split_condition.pdf"),
       split_plot, width = 12, height = 6)

saveRDS(obj, file.path(OUTPUT_DIR, "objects", "04_clustered.rds"))

message("=== Step 4 complete ===")
message(ncol(obj), " cells in ", length(levels(obj$seurat_clusters)), " clusters")
message("Review 04_elbow_plot.pdf; if N_DIMS looks wrong, edit it and re-run.")
message("Next: Rscript 05_cell_annotation.R ", DATASET)
