# Preprocessing notes and known gotchas

Practical notes on the Cell Ranger step. Read this if something does not behave
the way the README implies it should.

## Things to verify before submitting jobs

1. **Hashtag to sample assignments.** This is the single most common source of
   silently wrong results. If two hashtags are swapped in `config/samples.tsv`,
   every downstream comparison is wrong and nothing in the pipeline will warn
   you. Check the assignments against the original submission paperwork, not
   against a summary someone typed up later. Transcribing from merged cells in
   a Word or Excel table is especially error-prone.

2. **Excluded specimens.** If a specimen was dropped before sequencing, leave
   its row out of `samples.tsv` entirely. Its hashtag will simply go unused and
   Cell Ranger will not report it. Pools can legitimately have different sample
   counts.

3. **SLURM account and QoS.** Set `SLURM_ACCOUNT` and `SLURM_QOS` in
   `config/config.env` to the allocation the project should bill to. Jobs
   submitted to the wrong account may be rejected or charged to someone else.

4. **`CREATE_BAM`.** Defaults to `false` in this template. BAMs are large and
   most downstream work does not need them. Set it to `true` only if you plan
   to do variant calling or inspect alignments in IGV.

5. **Feature reference sequences.** `refs/feature_ref.csv` ships with the first
   four TotalSeq-C hashtag sequences as an example. If your panel differs and
   you do not update this file, you will get zero feature counts and
   demultiplexing will fail. Take the sequences from the vendor datasheet.

## Library types

`LIBRARY_TYPES` in `config/config.env` drives which sections appear in the
generated `multi_config_<POOL>.csv`:

| Library | `[libraries]` feature_types | Also adds |
|---|---|---|
| `GEX` | `Gene Expression` | required, always present |
| `CSP` | `Antibody Capture` | `[feature]` section |
| `TCR` | `VDJ-T` | `[vdj]` section |
| `BCR` | `VDJ-B` | `[vdj]` section |

The `[vdj]` section is written once if either TCR or BCR is present. Removing a
library type from the list removes its line from the config and its symlink
directory; no other edits are needed.

## CSP libraries: hashtags, antibodies, or both

A CSP (cell surface protein) library can carry:

- **hashtags only** — used purely to demultiplex pooled specimens
- **surface antibodies only** — quantified per cell, no demultiplexing
- **both**

All of them are `Antibody Capture` in the feature reference. The difference is
whether the ids appear in the `[samples]` section of the multi config, which
this template controls with `USE_HASHTAG_DEMUX`.

If you are not sure which you have, tally the feature barcodes present in the
CSP R2 reads and compare against your vendor's barcode list. A hashtag-only
library shows a small number of barcodes (one per hashtag) dominating the
counts.

Set `USE_HASHTAG_DEMUX=false` when the CSP library has no hashtags. The
`[samples]` section is then omitted and Cell Ranger produces one output per
pool rather than one per specimen.

## GEX spanning multiple flowcells or S-numbers

It is common for the gene expression library to be sequenced across more lanes
and flowcells than the CSP and V(D)J libraries, since GEX needs far more reads
per cell. That means one flat `gex_fastqs` symlink directory can contain files
with different S-numbers for the same pool.

Recent Cell Ranger versions accept this. If your version objects, list the
per-flowcell directories as a comma-separated `fastqs` value in the
`[libraries]` section instead of pointing at a single directory.

## Symlink counts do not match expectations

`01_make_symlinks.bash` matches files by the pattern
`*-<TAG>_*_R[12]_001.fastq.gz`. Sequencing cores do not all use the same naming
convention. If a library links zero files, print a few real filenames:

```bash
find "$RAW_FASTQ_DIR" -name '*.fastq.gz' | head
```

and adjust the `find` pattern in that script to match. Expect `2 x lanes x pools`
files per library type.

## Disk space

Cell Ranger output is large: budget roughly 20-50 GB per pool without BAMs, and
several times that with `CREATE_BAM=true`. The references add another 15-30 GB.
Check your quota before submitting a large batch of pools, since a job that
runs out of space partway through leaves an incomplete run directory that has
to be deleted and started again.

## Reruns

`04_submit_all.bash` skips any pool that already has an `outs` directory, so it
is safe to rerun after a partial failure. To force a pool to rerun, delete its
`<POOL>_run/<POOL>_multi/` directory first.

Cell Ranger also refuses to start if a run directory exists but is incomplete.
If a job was killed mid-run, delete the run directory rather than trying to
resume it.

## Suggested next steps not covered by this template

The R pipeline here takes you through QC, clustering, annotation, per-cluster
differential expression and integration across groups. Analyses deliberately
left out, because they are too study-specific to template usefully:

- **Subset re-clustering.** Pulling out one cell type and re-clustering it at
  higher resolution to resolve subpopulations.
- **Pseudobulk differential expression.** Aggregating counts per sample and
  running DESeq2 or edgeR. This is generally more statistically defensible than
  per-cell tests when you have real biological replicates, because per-cell
  tests treat cells as independent replicates and badly inflate significance.
- **Clonotype analysis.** The V(D)J contigs are loaded and attached to the
  Seurat object, but clonal expansion, diversity and overlap analyses are not
  implemented. Consider `scRepertoire` for this.
- **Trajectory and cell-cell communication analyses.**
