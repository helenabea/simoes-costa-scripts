# Core Infrastructure – Simoes-Costa Lab

This repository contains core analysis scripts and computational utilities used in the Simoes-Costa Lab for the processing and analysis of next-generation sequencing (NGS) data.

The goal of this repository is to provide **transparent, easy-to-understand scripts** that can be used and adapted by lab members for common genomics workflows.

---

## Available Pipelines

The repository currently provides scripts supporting the following workflows:

### RNA-seq
RNA-seq alignment and processing pipelines using:

- **HISAT2**
- **STAR**

Both paired-end and single-end configurations are supported.

### ATAC-seq
Basic processing workflow for ATAC-seq experiments including alignment and peak calling.

### CUT&RUN
Scripts for CUT&RUN data processing and alignment.

---

## Usage

1. Edit the **configuration section** at the top of the relevant script.
2. Specify genome references and FASTQ inputs.
3. Run the script from the command line.

---

## Notes

Each script contains a **configuration section** where users can specify:

- reference genome
- annotation files
- input FASTQ files
- output directories

---

## Requirements

These workflows assume a Linux environment with common bioinformatics tools installed.

Some dependencies include:

- fastqc
- bowtie2
- HISAT2
- STAR
- samtools
- macs2
- featureCounts

---

## Maintainer

Computational organization and infrastructure currently maintained by:

Helena B. Conceição  
Postdoctoral Researcher – Computational Genomics
