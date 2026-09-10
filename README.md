# Core Pipelines – Simoes-Costa Lab

This repository contains core analysis scripts and computational utilities used in the Simoes-Costa Lab for the processing and analysis of next-generation sequencing (NGS) data.

The goal of this repository is to provide **transparent, easy-to-understand scripts** that can be used and adapted by lab members for common genomics workflows.

---

## Available Pipelines

The repository currently provides scripts supporting the following workflows:

- **RNA-seq** alignment and processing pipelines using:
  - **HISAT2**
  - **STAR**
- **ATAC-seq** alignment and processing pipelines including peak calling.
- **CUT&RUN** alignment and processing pipelines including peak calling.

Both paired-end and single-end configurations are supported for **RNA-seq**.

---

## Choose a pipeline

| Assay | Read layout | Script |
|---|---|---|
| RNA-seq | Paired-end | `STAR_rnaseq_PE.sh` or `HiSat2_rnaseq_PE.sh` |
| RNA-seq | Single-end | `STAR_rnaseq_SE.sh` or `HiSat2_rnaseq_SE.sh` |
| ATAC-seq | Paired-end | `ATAC_pipeline.sh` |
| CUT&RUN histone marks or transcription factors | Paired-end | `CutRun_pipeline.sh` |

## Usage

1. Copy the selected script into a new, empty analysis directory.
2. Edit only the section marked **USER CONFIGURATION** at the top.
3. Open a new tmux/screen session.
4. Run the pipeline with `bash SCRIPT_NAME.sh`.

Each pipeline creates its output folders in the current working directory. Using a
new directory for every analysis prevents results from different runs from being
mixed together or overwritten.

## Parameters users should edit

- `GENOME`: one of `galgal7`, `galgal6`, `hg38`, or `mm39`.
- `THREADS`: number of CPU threads available to the job. Use less than 10.
- `ORIGINAL_FASTQ_PATH`: absolute path to the directory containing the original
  compressed FASTQ files.
- `SAMPLE_ID`: one or more unique prefixes from the sample sheet. Quote every
  entry, for example `SAMPLE_ID=("sample_rep1" "sample_rep2")`.
- `R1_PATTERN` and `R2_PATTERN`: exact filename suffixes used to identify reads.
  Paired-end scripts use both; single-end scripts use only `R1_PATTERN`. The text
  before the suffix becomes the output sample name, so it must be unique.
- `MAPQ_THRESHOLD` (ATAC-seq and CUT&RUN): minimum alignment mapping quality
  retained for downstream analysis. The defaults are `30` for ATAC-seq and `0`
  for CUT&RUN. CUT&RUN of Histone marks would benefit from setting MAPQ threshold to `20`.
- `REMOVE_DUPLICATES` (CUT&RUN): read/PCR duplicate removal is now optional.
  Duplicates are always marked and measured first.
- `TRACK_BIN_SIZE`: BigWig resolution in base pairs. RNA-seq pipelines produce
  one combined-strand CPM track per sample. ATAC-seq and CUT&RUN produce both
  `.CPM.bw` and `.RPGC.bw` tracks in the `BW` directory. RPGC uses the
  reference-specific `GENOME_SIZE` defined in the genome-settings section.
- `FEATURECOUNTS_STRANDED` (RNA-seq): `0` for unstranded, `1` for forward, or
  `2` for reverse-stranded libraries. Confirm this from the library preparation
  method rather than guessing.
- `FEATURECOUNTS_MODE` (RNA-seq): use `"matrix"` to create one count table
  containing every processed sample, or `"per_sample"` to create a separate
  count file for each sample. The default is `"matrix"`.
- `HISAT_STRANDNESS` (HISAT2 only): `F`/`R` for single-end or `FR`/`RF` for
  paired-end libraries.
- `KEEP_TRIMMED_FASTQ`: set to `true` to retain trimmed reads or `false` to
  remove them after successful downstream processing.
- `START_AT`: controls full and resumed runs. Keep `"all"` to run the
  complete pipeline. RNA-seq pipelines also accept `"trim"`, `"align"`, or
  `"count"`. ATAC-seq and CUT&RUN accept `"trim"`, `"align"`, `"duplicates"`,
  `"tracks"`, or `"peaks"`. The selected stage and every downstream stage run.

For example, to repeat alignment and counting after an interrupted RNA-seq run:

```bash
START_AT="align"
```

Resume the pipeline from the same analysis directory so it can find outputs from
the earlier stages. Starting at `align` requires trimmed FASTQ files; therefore,
set `KEEP_TRIMMED_FASTQ=true` on long-running jobs where a later restart may be
needed.

Genome index and annotation paths are maintained in the protected genome-settings
section of each script. If the shared reference directories move, update those
paths consistently across all pipeline scripts.

---

## Maintainer

Computational organization and infrastructure currently maintained by:

Helena B. Conceição
Postdoctoral Researcher – Computational Genomics
