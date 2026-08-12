# =============================================================================
# Step 7: Copy results to the web publishing directory
# =============================================================================
# Gathers the plots, tables and Cell Ranger web summaries for one group into
# WEB_PUBLISH_DIR so collaborators can browse them.
#
# Skipped if WEB_PUBLISH_DIR is still the placeholder value in config.env,
# so the pipeline works fine if you do not publish results to the web.
#
# Run:
#     Rscript 07_prepare_report.R [GROUP_NAME]
# =============================================================================

# Find and load the shared config by walking up to the repo root.
.p <- normalizePath(getwd(), mustWork = FALSE)
while (!file.exists(file.path(.p, "02_R", "R", "config.R")) && dirname(.p) != .p) {
  .p <- dirname(.p)
}
if (!file.exists(file.path(.p, "02_R", "R", "config.R"))) {
  stop("Could not find 02_R/R/config.R. Run this script from inside the repository.")
}
source(file.path(.p, "02_R", "R", "config.R"))

message("=== Step 7: Preparing report files for group '", DATASET, "' ===")

if (grepl("^/path/to/", CFG$WEB_PUBLISH_DIR %||% "/path/to/")) {
  message("WEB_PUBLISH_DIR is not configured in config/config.env; skipping.")
  message("Results remain in: ", OUTPUT_DIR)
  quit(save = "no", status = 0)
}

for (sub in c("plots", "tables", "web_summaries")) {
  dir.create(file.path(WEB_DIR, sub), recursive = TRUE, showWarnings = FALSE)
}

# ---- analysis outputs ----
n_copied <- 0
for (sub in c("plots", "tables")) {
  src <- file.path(OUTPUT_DIR, sub)
  if (!dir.exists(src)) next
  files <- list.files(src, full.names = TRUE)
  if (length(files) == 0) next
  ok <- file.copy(files, file.path(WEB_DIR, sub), overwrite = TRUE)
  n_copied <- n_copied + sum(ok)
  message("  ", sub, ": copied ", sum(ok), " files")
}

# ---- Cell Ranger web summaries ----
# These are the per-sample QC pages; useful to publish alongside the analysis
# so collaborators can check the upstream data quality themselves.
per_sample_dir <- file.path(CELLRANGER_OUT, "per_sample_outs")
if (dir.exists(per_sample_dir)) {
  sample_dirs <- list.dirs(per_sample_dir, full.names = TRUE, recursive = FALSE)
  for (sd in sample_dirs) {
    src <- file.path(sd, "web_summary.html")
    if (file.exists(src)) {
      file.copy(src,
                file.path(WEB_DIR, "web_summaries",
                          paste0(basename(sd), "_web_summary.html")),
                overwrite = TRUE)
      n_copied <- n_copied + 1
    }
  }
  message("  web_summaries: copied ", length(sample_dirs), " summaries")
} else {
  message("  web_summaries: Cell Ranger output not found, skipping")
}

message("=== Step 7 complete ===")
message(n_copied, " files copied to: ", WEB_DIR)
if (nzchar(CFG$WEB_PUBLISH_URL %||% "")) {
  message("Should be visible at: ", WEB_URL)
}
message("Next: render the report with  Rscript -e \"rmarkdown::render('08_report.Rmd')\"")
