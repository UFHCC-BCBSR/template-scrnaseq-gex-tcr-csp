# =============================================================================
# Step 5: Cluster markers and cell type annotation aids
# =============================================================================
# Finds the genes that distinguish each cluster, and plots canonical cell type
# markers so you can work out what each cluster is.
#
# This step does NOT assign cell type labels automatically. Automated labelling
# is unreliable without a well-matched reference, and a wrong label propagates
# silently into every downstream result. Use the outputs here to annotate
# manually; see the EDIT block at the bottom.
#
# Run:
#     Rscript 05_cell_annotation.R [GROUP_NAME]
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

obj <- readRDS(file.path(OUTPUT_DIR, "objects", "04_clustered.rds"))

message("=== Step 5: Cell annotation for group '", DATASET, "' ===")

# ---- cluster marker genes ----
message("Finding cluster markers (this takes a few minutes)...")
cluster_markers <- FindAllMarkers(
  obj,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)
write_csv(cluster_markers, file.path(OUTPUT_DIR, "tables", "05_cluster_markers.csv"))

top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write_csv(top_markers, file.path(OUTPUT_DIR, "tables", "05_top_markers_per_cluster.csv"))

message("Top markers per cluster:")
print(top_markers %>% select(cluster, gene, avg_log2FC, p_val_adj), n = 100)

# ---- canonical marker plots ----
# CELL_TYPE_MARKERS comes from config.R and is species-appropriate. Add your
# own tissue-specific markers there, or override the list here.
message("\nPlotting canonical cell type markers...")

for (cell_type in names(CELL_TYPE_MARKERS)) {
  markers   <- CELL_TYPE_MARKERS[[cell_type]]
  available <- markers[markers %in% rownames(obj)]
  missing   <- setdiff(markers, available)

  if (length(missing) > 0) {
    message("  ", cell_type, ": not in data: ", paste(missing, collapse = ", "))
  }
  if (length(available) == 0) next

  p <- FeaturePlot(obj, features = available, ncol = 2) &
    theme(plot.title = element_text(size = 10))
  p <- p + plot_annotation(
    title = cell_type,
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
  )
  ggsave(file.path(OUTPUT_DIR, "plots", paste0("05_markers_", cell_type, ".pdf")),
         p, width = 12, height = 8)
}

# A dot plot of all markers at once is often the fastest way to read off
# which cluster is which.
all_markers <- unique(unlist(CELL_TYPE_MARKERS))
all_markers <- all_markers[all_markers %in% rownames(obj)]
if (length(all_markers) > 0) {
  dp <- DotPlot(obj, features = all_markers) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
    ggtitle("Canonical markers by cluster")
  ggsave(file.path(OUTPUT_DIR, "plots", "05_marker_dotplot.pdf"),
         dp, width = max(12, length(all_markers) * 0.25), height = 8)
}


# ---- EDIT: manual cluster annotation ----------------------------------------
# Once you have reviewed the marker plots and the dot plot, map each cluster
# number to a cell type here, then re-run this script.
#
# Uncomment and fill in. Every cluster present in the data must have an entry.
#
# cluster_annotations <- c(
#   "0" = "CD4 T cells",
#   "1" = "CD8 T cells",
#   "2" = "B cells",
#   "3" = "Macrophages"
# )
#
# obj$cell_type <- unname(cluster_annotations[as.character(obj$seurat_clusters)])
# if (any(is.na(obj$cell_type))) {
#   stop("These clusters have no annotation: ",
#        paste(sort(unique(obj$seurat_clusters[is.na(obj$cell_type)])), collapse = ", "))
# }
#
# p <- DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE) +
#   ggtitle("Annotated cell types")
# ggsave(file.path(OUTPUT_DIR, "plots", "05_umap_annotated.pdf"), p, width = 12, height = 8)

if (!"cell_type" %in% colnames(obj@meta.data)) {
  message("\nNOTE: clusters are not yet annotated. Downstream steps will use ",
          "cluster numbers.\n      Fill in the 'EDIT: manual cluster annotation' ",
          "block in this script to add cell type labels.")
}

saveRDS(obj,             file.path(OUTPUT_DIR, "objects", "05_annotated.rds"))
saveRDS(cluster_markers, file.path(OUTPUT_DIR, "objects", "05_cluster_markers.rds"))

message("=== Step 5 complete ===")
message("Next: Rscript 06_de.R ", DATASET)
