# scRNA-seq analysis template: GEX + TCR + CSP

A starting point for 10x Genomics single-cell RNA-seq projects that combine
gene expression with V(D)J and cell surface protein / hashtag libraries.

It covers the whole path from raw fastq files to an HTML report: read quality
checks, Cell Ranger, then a Seurat pipeline for QC, clustering, annotation,
differential expression, and integration across groups.

This README is written to be followed step by step by someone who has not used
the command line before. Read it top to bottom the first time. Every command
you need to type is shown in a grey box. Type or paste one box at a time and
press Enter.

---

## Using this template

Click **Use this template** at the top of the GitHub page to create your own
copy, then clone it onto your cluster:

```bash
git clone https://github.com/YOUR-ORG/YOUR-NEW-REPO.git
cd YOUR-NEW-REPO
```

Everything project-specific lives in two files in `config/`. You should not
need to edit the pipeline scripts to run a standard analysis.

### What this template assumes

- 10x Genomics data processed with `cellranger multi`
- One or more sequencing pools, optionally multiplexed with hashtag antibodies
- A Linux cluster running the SLURM scheduler
- Mouse or human samples (other species need new reference URLs)

If your setup differs, the pipeline is still a reasonable starting point, but
expect to adjust the preprocessing scripts.

---

## Part 1: Background

### Library types

A typical experiment of this kind sequences several library types from the same
cells, each answering a different question:

| Library | What it measures | Typical depth |
|---|---|---|
| GEX | gene expression | ~30,000 reads per cell |
| TCR | T cell receptor sequences | ~5,000 reads per cell |
| BCR | B cell receptor sequences | ~5,000 reads per cell |
| CSP | cell surface proteins and/or hashtags | ~5,000 reads per cell |

The default configuration is GEX + TCR + CSP. Add or remove types by editing
`LIBRARY_TYPES` in `config/config.env`; the scripts adapt automatically.

### Pooling and hashtags

Specimens are often pooled before sequencing to save money and reduce batch
effects, then separated again computationally. Hashtag antibodies (for example
TotalSeq-C) tag each specimen's cells with a unique barcode that is read out in
the CSP library. Cell Ranger uses these to split a pool back into individual
specimens.

If you did not pool, set `USE_HASHTAG_DEMUX=false` and give each specimen its
own pool in the sample sheet.

### What preprocessing produces

The goal is to turn raw sequencing files into per-cell gene, receptor and
protein counts, and to split each pool back into its individual specimens. The
tool that does this is Cell Ranger, run once per pool.

---

## Part 2: Getting onto your cluster

> **Adapt this section to your institution.** The commands below are generic;
> the hostnames and web addresses are not. Ask your local research computing
> group for the equivalents, or check their documentation.

Most university clusters offer two ways in.

### Option A: a web portal (easiest for beginners)

Many clusters run [Open OnDemand](https://openondemand.org/), which gives you a
file browser and a terminal inside your web browser, and can display the HTML
quality reports this pipeline produces. Look for a URL like
`https://ondemand.your-institution.edu`.

Once logged in:

- **Terminal**: Clusters menu, then Shell Access.
- **Files**: Files menu, then Home Directory.

### Option B: ssh from your own machine

On a Mac, open the Terminal app. On Windows, open Windows Terminal or
PowerShell. Then:

```bash
ssh your_username@login.your-cluster.edu
```

Type your password when asked. Nothing appears on screen as you type the
password; that is normal. Press Enter.

---

## Part 3: Command line basics

The command line is a way to tell the computer what to do by typing instead of
clicking. These are the only commands this guide needs.

Show which folder you are currently in:

```bash
pwd
```

List the files and folders here:

```bash
ls
```

Move into a folder:

```bash
cd 01_preprocessing
```

Move up one folder:

```bash
cd ..
```

Print a small text file to the screen:

```bash
cat config/config.env
```

View a longer file one screen at a time (space to page down, `q` to quit):

```bash
less README.md
```

Tips:

- After typing part of a name, press Tab and the computer will finish it. This
  saves typing and avoids spelling mistakes.
- Press the up arrow to bring back the previous command.
- Commands and file names are case sensitive. `Ls` is not the same as `ls`.

---

## Part 4: Configure your project

This is the only part you must edit by hand, and everything downstream depends
on getting it right.

### Step 4.1: Edit `config/config.env`

```bash
nano config/config.env
```

In nano, use the arrow keys to move, type to edit, press Control-O then Enter
to save, and Control-X to exit.

At minimum, set:

| Setting | What to put |
|---|---|
| `PROJECT_ROOT` | absolute path to this checkout on the cluster |
| `RAW_FASTQ_DIR` | where your sequencing core delivered the fastq files |
| `SPECIES` | `mouse` or `human` |
| `LIBRARY_TYPES` | which libraries you sequenced |
| `SLURM_ACCOUNT` / `SLURM_QOS` | the allocation to bill jobs to |

The file documents every setting inline. The scripts refuse to run while the
placeholder values are still in place, so you cannot forget by accident.

### Step 4.2: Edit `config/samples.tsv`

This describes your experiment: one row per final demultiplexed specimen.

```bash
nano config/samples.tsv
```

The shipped file contains a worked example with two pools and four hashtags
each. Replace it entirely with your own samples. The columns are documented in
comments at the top of the file.

> **This is the most important file in the project.** If two hashtags are
> swapped here, every downstream result will be wrong and nothing will warn
> you. Check the assignments against the original submission paperwork, not
> against a summary someone typed up later.

Check that it reads correctly before continuing:

```bash
column -t -s $'\t' config/samples.tsv | grep -v '^#'
```

---

## Part 5: Check the raw data quality first

Before running Cell Ranger, confirm the sequencing itself worked.

### Step 5.1: Build the list of fastq files

```bash
cd 01_preprocessing/01_fastqc
bash 00_make_sample_list.bash
```

This searches your delivery directory and writes `sample_fastqs.tsv`. It prints
how many read pairs it found per library type. **Compare these against your
sequencing core's delivery manifest.** A library sequenced across four lanes
should show four pairs per pool.

### Step 5.2: Run FastQC

```bash
bash 02_submit_fastqc.bash
```

This submits one job per fastq pair to the cluster. Watch progress with:

```bash
squeue -u $USER
```

When the jobs disappear from that list, they have finished.

### Step 5.3: Summarise with MultiQC

```bash
bash 03_multiqc.bash
```

Then open `reports/multiqc_report.html` in a browser. Using a web portal: go to
the Files menu, navigate to the folder, click the three dots next to the file
and choose View.

### Step 5.4: What to check, and what is a red flag

Single cell data looks different from ordinary sequencing, so FastQC marks
several things orange or red that are completely normal here. Do not be alarmed
by colour alone.

**Expected, and NOT a problem for 10x data:**

- *Per Base Sequence Content* flagged at the start of the read. The cell
  barcode and primer regions are not random, so this is always flagged.
- *Sequence Duplication Levels* flagged high. Single cell libraries are
  amplified and targeted, so high duplication is expected.
- *Overrepresented Sequences*. Poly-A stretches and, for CSP libraries, the
  hashtag barcodes themselves show up here.
- Read 1 being short, around 26 to 28 bases. Read 1 holds the cell barcode and
  UMI, not gene sequence.

**Genuine red flags, worth stopping for:**

- *Per Base Sequence Quality* dropping into the red across most of Read 2. A
  little decline at the very end of a read is fine; a whole read that is low
  quality is not.
- *Per Sequence Quality Scores* centred low. Most reads should have a mean
  quality above 30.
- One sample with far fewer reads than others of the same library type. Compare
  the read counts in the General Statistics table at the top. A pool with a
  tiny fraction of the others' reads may be a failed library.
- A sample missing entirely from a report where you expected it.
- *Adapter Content* climbing to a large fraction of the read. Cell Ranger trims
  some adapter, so a small amount is fine; a read that is mostly adapter is not.

Write down anything in the second list and discuss it before proceeding.

---

## Part 6: Run Cell Ranger

```bash
cd ../02_cellranger
```

### Step 6.1: Build the references (once, ever)

Downloads the genome and V(D)J references for your species. This is a large
download and can take a while.

```bash
bash 02_make_refs.bash
```

The script skips anything already present, so it is safe to re-run. If your
cluster already provides shared 10x references, you can skip this and point at
those instead.

### Step 6.2: Gather the fastq files

Sequencing cores deliver fastqs in deep nested folders; Cell Ranger wants them
flat. This creates one tidy directory of symlinks per library type without
copying the large files.

```bash
bash 01_make_symlinks.bash
```

It prints how many files it linked per library. **Check these numbers.** You
should see `2 x lanes x pools` files per library type. If they look wrong, stop
and investigate before continuing.

### Step 6.3: Generate the run configurations

```bash
bash 03_make_per_pool_runs.bash
```

This writes one folder per pool, each with a Cell Ranger configuration file and
a job script. Look at one before continuing:

```bash
cat POOL1_run/multi_config_POOL1.csv
```

You should see your library types listed, and a `[samples]` section listing
each specimen against its hashtag. Check the specimen names and hashtags match
your sample sheet.

### Step 6.4: Submit the jobs

```bash
bash 04_submit_all.bash
```

Each submission prints a job ID. A Cell Ranger run typically takes several
hours. Watch them with:

```bash
squeue -u $USER
```

To follow a job's log as it runs (Control-C to stop watching):

```bash
tail -f POOL1_run/logs/cellranger_POOL1.*.out
```

### Step 6.5: Check the results

```bash
bash 05_check_runs.bash
```

This reports which runs finished and lists the web summaries to open. For each
one, check that:

- the estimated number of cells is plausible and not near zero
- the multiplexing section shows the expected number of specimens for that pool
- there are no large red warning banners

If a run failed, the script prints the end of its error log.

---

## Part 7: Downstream analysis in R

```bash
cd ../../02_R/per_sample_analysis
```

### Running the pipeline

Run every step for every group:

```bash
bash run_all.bash
```

Or one group at a time:

```bash
bash run_all.bash TISSUE_A
```

Or step by step, which is what you want while you are still tuning parameters:

```bash
cd R
Rscript 01_load_data.R TISSUE_A
Rscript 02_qc_filtering.R TISSUE_A
```

### The steps

| Script | What it does |
|---|---|
| `01_load_data.R` | Reads Cell Ranger output into a Seurat object |
| `02_qc_filtering.R` | Computes QC metrics and removes low-quality cells |
| `03_vdj_integration.R` | Attaches TCR/BCR data to the object |
| `04_norm_and_cluster.R` | Normalisation, PCA, clustering, UMAP |
| `05_cell_annotation.R` | Cluster markers and cell type marker plots |
| `06_de.R` | Differential expression between conditions |
| `07_prepare_report.R` | Copies results to the web directory |
| `08_report.Rmd` | Renders the HTML report |
| `09_loupe.R` | Exports a `.cloupe` file for Loupe Browser |

### Parameters you should review

The defaults are reasonable starting points, not answers. Three in particular
deserve your attention:

- **QC thresholds** at the top of `02_qc_filtering.R`. Look at the "before
  filtering" plots first, then set cutoffs that match what you see.
- **`N_DIMS`** in `04_norm_and_cluster.R`. Check `04_elbow_plot.pdf` and use
  the number of principal components where the curve flattens.
- **Cluster annotation** in `05_cell_annotation.R`. There is an `EDIT` block
  near the bottom for mapping cluster numbers to cell types. Automated
  labelling is deliberately not used, because a wrong label propagates silently
  into every downstream result.

### Combining groups

Once each group has been through the per-sample pipeline:

```bash
cd ../integrated_analysis/R
Rscript 01_integrate_data.R
Rscript 02_loupe.R
```

This merges the groups and applies Harmony batch correction. Check
`integration_by_group.pdf`: the groups should overlap. If they form separate
islands, the correction has not worked and clusters may reflect batch rather
than biology.

---

## Part 8: Save your work with git

git keeps a history of the project and syncs it to GitHub so nothing is lost.

```bash
git add config/ 
git commit -m "Configure project for <your experiment>"
git push
```

Large outputs (Cell Ranger results, references, symlinks, `.rds` objects,
`.cloupe` files) are excluded by `.gitignore` on purpose. Only the small
configuration and script files belong in git.

---

## Repository structure

```
config/
  config.env              paths, species, library types, cluster settings
  samples.tsv             one row per specimen
  load_config.bash        shared loader for the shell scripts

01_preprocessing/
  01_fastqc/              read quality checks
  02_cellranger/          references, symlinks, per-pool runs

02_R/
  R/config.R              shared config loader for the R scripts
  per_sample_analysis/    the 9-step pipeline, run once per group
  integrated_analysis/    combine groups with Harmony

03_reports/               progress update slide template

docs/
  preprocessing_notes.md  gotchas, library type details, next steps
```

---

## Troubleshooting

**A script says a config value is still a placeholder.** Edit
`config/config.env` and replace the `/path/to/...` values.

**`01_make_symlinks.bash` links zero files.** Your core uses a different naming
convention. Print some real filenames with
`find "$RAW_FASTQ_DIR" -name '*.fastq.gz' | head` and adjust the pattern in the
script.

**Cell Ranger recovers fewer specimens than expected.** Check the hashtag
assignments in `config/samples.tsv` and the barcode sequences in
`01_preprocessing/02_cellranger/refs/feature_ref.csv`.

**`percent.mt` is zero for every cell.** `SPECIES` does not match your
reference. Mouse mitochondrial genes are `mt-*`, human are `MT-*`.

**An R script cannot find `config.R`.** Run it from inside the repository; the
scripts locate the config by walking up the directory tree.

More detail in `docs/preprocessing_notes.md`.

---

## License

MIT. See [LICENSE](LICENSE).
