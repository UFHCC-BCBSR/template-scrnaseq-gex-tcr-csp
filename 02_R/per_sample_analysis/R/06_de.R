# =============================================================================
# Step 6: Differential expression between conditions
# =============================================================================
# Compares conditions overall and within each cluster.
#
# STATISTICAL CAVEAT: these are per-cell tests. They treat every cell as an
# independent observation, which it is not - cells from one animal or donor are
# correlated. This inflates significance, sometimes dramatically. Using
# sample_id as a latent variable (as below) mitigates but does not remove the
# problem.
#
# If you have three or more biological replicates per condition, a pseudobulk
# analysis (aggregate counts per sample, then DESeq2 or edgeR) is more
# defensible and should be preferred for anything going into a paper. See
# docs/preprocessing_notes.md.
#
# Run:
#     Rscript 06_de.R [GROUP_NAME]
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
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


# ---- parameters -------------------------------------------------------------
LOGFC_THRESHOLD <- 0.25
MIN_CELLS_PER_GROUP <- 10   # skip a cluster if either condition has fewer

# MAST accounts for the cellular detection rate and allows latent variables.
# Falls back to the Wilcoxon test if MAST is not installed.
DE_TEST <- if (requireNamespace("MAST", quietly = TRUE)) "MAST" else "wilcox"


obj <- readRDS(file.path(OUTPUT_DIR, "objects", "05_annotated.rds"))

message("=== Step 6: Differential expression for group '", DATASET, "' ===")
message("Test: ", DE_TEST)
if (DE_TEST == "wilcox") {
  message("NOTE: MAST is not installed, using the Wilcoxon test instead. ",
          "Install MAST for detection-rate-aware testing.")
}

# ---- work out what to compare ----
# Conditions come from config/samples.tsv rather than being hard-coded.
conditions <- unique(obj$condition)
conditions <- conditions[!is.na(conditions)]

if (length(conditions) < 2) {
  stop("Group '", DATASET, "' has fewer than two conditions (",
       paste(conditions, collapse = ", "), "), so there is nothing to compare.\n",
       "Check the 'condition' column in config/samples.tsv.")
}

# All pairwise comparisons. With exactly two conditions this is a single
# comparison, which is the common case.
comparisons <- combn(sort(conditions), 2, simplify = FALSE)
message("Comparisons: ",
        paste(sapply(comparisons, function(x) paste(x[2], "vs", x[1])), collapse = "; "))

# Use cell type labels if step 5 produced them, otherwise cluster numbers.
cluster_var <- if ("cell_type" %in% colnames(obj@meta.data)) "cell_type" else "seurat_clusters"
message("Grouping cells by: ", cluster_var)

latent <- if (DE_TEST == "MAST" && length(unique(obj$sample_id)) > 1) "sample_id" else NULL


#' Run one differential expression comparison, returning NULL if not possible
run_de <- function(object, ident_1, ident_2, label) {
  n1 <- sum(object$condition == ident_1, na.rm = TRUE)
  n2 <- sum(object$condition == ident_2, na.rm = TRUE)

  if (n1 < MIN_CELLS_PER_GROUP || n2 < MIN_CELLS_PER_GROUP) {
    message("  ", label, ": skipped (", ident_1, " n=", n1, ", ", ident_2, " n=", n2, ")")
    return(NULL)
  }

  Idents(object) <- "condition"
  res <- tryCatch(
    FindMarkers(object,
                ident.1         = ident_1,
                ident.2         = ident_2,
                test.use        = DE_TEST,
                latent.vars     = latent,
                logfc.threshold = LOGFC_THRESHOLD),
    error = function(e) {
      message("  ", label, ": failed (", conditionMessage(e), ")")
      NULL
    }
  )
  if (is.null(res) || nrow(res) == 0) return(NULL)

  res$gene <- rownames(res)
  n_sig <- sum(res$p_val_adj < 0.05, na.rm = TRUE)
  message("  ", label, ": ", nrow(res), " genes tested, ", n_sig, " with adj p < 0.05")
  res
}


# ---- overall comparisons ----
message("\nOverall (all cells):")
for (cmp in comparisons) {
  ref <- cmp[1]; alt <- cmp[2]
  tag <- paste0(alt, "_vs_", ref)

  res <- run_de(obj, alt, ref, tag)
  if (!is.null(res)) {
    write_csv(res, file.path(OUTPUT_DIR, "tables",
                             paste0("06_de_overall_", tag, ".csv")))

    # Volcano plot
    res$significant <- res$p_val_adj < 0.05 & abs(res$avg_log2FC) > 0.5
    v <- ggplot(res, aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300),
                         colour = significant)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
      labs(title = paste("Differential expression:", alt, "vs", ref),
           subtitle = paste("Group:", DATASET),
           x = "log2 fold change", y = "-log10 adjusted p") +
      theme_bw() + theme(legend.position = "none")
    ggsave(file.path(OUTPUT_DIR, "plots", paste0("06_volcano_", tag, ".pdf")),
           v, width = 8, height = 6)
  }
}

# ---- per-cluster comparisons ----
message("\nWithin each ", cluster_var, ":")
for (cmp in comparisons) {
  ref <- cmp[1]; alt <- cmp[2]

  for (cl in sort(unique(obj@meta.data[[cluster_var]]))) {
    cells <- rownames(obj@meta.data)[obj@meta.data[[cluster_var]] == cl]
    if (length(cells) < 2 * MIN_CELLS_PER_GROUP) next

    subset_obj <- subset(obj, cells = cells)
    tag <- paste0(alt, "_vs_", ref, "_", gsub("[^A-Za-z0-9]+", "_", as.character(cl)))

    res <- run_de(subset_obj, alt, ref, paste0("  ", cluster_var, " ", cl))
    if (!is.null(res)) {
      res$cluster <- cl
      write_csv(res, file.path(OUTPUT_DIR, "tables", paste0("06_de_", tag, ".csv")))
    }
  }
}

# ---- composition summary ----
# Shifts in how many cells fall in each cluster are often as biologically
# interesting as changes in expression within a cluster.
summary_stats <- obj@meta.data %>%
  group_by(condition, .data[[cluster_var]]) %>%
  summarise(
    n_cells      = n(),
    bcr_positive = sum(has_bcr, na.rm = TRUE),
    tcr_positive = sum(has_tcr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = condition,
              values_from = c(n_cells, bcr_positive, tcr_positive),
              names_sep = "_", values_fill = 0)
write_csv(summary_stats, file.path(OUTPUT_DIR, "tables", "06_cluster_composition.csv"))

comp_plot <- obj@meta.data %>%
  count(.data[[cluster_var]], condition) %>%
  ggplot(aes(x = factor(.data[[cluster_var]]), y = n, fill = condition)) +
  geom_col(position = "fill") +
  labs(title = "Cluster composition by condition", subtitle = paste("Group:", DATASET),
       x = cluster_var, y = "proportion of cells") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "plots", "06_cluster_composition.pdf"),
       comp_plot, width = 10, height = 6)

saveRDS(obj, file.path(OUTPUT_DIR, "objects", "06_final_analysis.rds"))

message("\n=== Step 6 complete ===")
message("Result tables: ", file.path(OUTPUT_DIR, "tables"))
message("Next: Rscript 07_prepare_report.R ", DATASET)
