# Microbial Signatures in Host-Unmapped Human RNA-seq Reads

## Overview

This project develops a reproducible bioinformatics workflow for investigating microbial signatures in host-unmapped reads from human RNA-seq datasets.

Human RNA-seq datasets contain reads that may fail to align to the human reference genome. Rather than discarding these host-unmapped reads, this project examines them further using taxonomic classification and abundance estimation to investigate potential microbial signatures.

The workflow is being developed in two stages. An initial single-sample pilot was used to establish and troubleshoot the local command-line workflow. The pipeline will subsequently be refined and scaled to analyse 15 human tissue RNA-seq samples from skin, lung and colon.

## Project Rationale

Human RNA-seq studies are primarily designed to investigate host gene expression. During alignment to a human reference genome, however, a proportion of reads may remain unmapped. These reads are often excluded from host-focused downstream analyses, but they may contain additional biological information that can be investigated further.

This project investigates whether host-unmapped RNA-seq reads can be characterise to identify potential microbial signatures in human tissue transcriptomic datasets. Taxonomic classification and abundance estimation are used to examine the microbial taxa represented within these reads while recognising that host-unmapped reads are not necessarily microbial in origin.

## Workflow Development

The workflow was first established using a single RNA-seq sample as a pilot. This initial implementation was used to test the major stages of the pipeline, identify computational and methodological challenges, and determine the refinements required before scaling the analysis to the full 15-sample dataset.

### Pilot Workflow

The pilot workflow processed a single-end RNA-seq sample through quality control, read preprocessing, host alignment, microbial taxonomic classification and abundance estimation. FastQC v0.12.1 was used for read quality assessment, fastp v1.3.7 for preprocessing, HISAT2 v2.2.3 for alignment against the human reference genome, Kraken2 v2.17.1 for taxonomic classification of host-unmapped reads, and Bracken v3.0.1 for species-level abundance estimation. MultiQC v1.35 was used to aggregate and summarise quality-control and pipeline reports.

Pilot workflow:

`Raw RNA-seq reads → FastQC → fastp → FastQC → HISAT2 → host-unmapped reads → Kraken2 → Bracken`

#### Pilot Sample

The pilot workflow was implemented using SRR7961304 (SKIN_3), a 50 bp single-end human skin RNA-seq sample included in the full 15-sample dataset. This sample was selected to establish and troubleshoot the local command-line implementation of the workflow before scaling the pipeline across the complete dataset. The pilot therefore serves as an initial implementation stage rather than a separate analysis from the full study.

### Pilot Observations and Limitations

The pilot successfully established the major stages of the local workflow but also identified areas requiring refinement before multi-sample analysis. In particular, Kraken2 classification using the Standard-Full database was computationally demanding on the available local hardware, resulting in a substantial runtime. The pilot therefore highlighted the need to evaluate a more efficient approach before scaling microbial classification to all 15 samples.

The pilot also informed several methodological refinements for the full analysis. This repository extends a microbial-signature analysis of host-unmapped human RNA-seq reads originally undertaken as part of an MSc research project into a reproducible local command-line workflow. The full 15-sample analysis will implement the original MSc preprocessing approach, using Falco for quality control and Trimmomatic for read preprocessing. Quality control will also be performed on the host-unmapped reads before they proceed to microbial taxonomic classification. In addition, the Bash pipeline will be expanded to support both single-end and paired-end RNA-seq data, applying the appropriate processing commands according to the sequencing layout of each sample.

### Refined 15-Sample Workflow

Building on the pilot, the refined workflow will process the complete 15-sample dataset comprising five skin, five lung and five colon tissue RNA-seq samples. The workflow will retain host alignment with HISAT2 followed by taxonomic classification and abundance estimation with Kraken2 and Bracken, while incorporating the preprocessing, quality-control and sequencing-layout refinements identified above.

Refined workflow:

`Raw RNA-seq reads → Falco → Trimmomatic → Falco → HISAT2 → host-unmapped reads → Falco → Kraken2 → Bracken`

## Dataset

The full analysis comprises 15 human tissue RNA-seq samples, with five biological samples representing each of three tissue groups: skin, lung and colon. The dataset includes both single-end and paired-end sequencing layouts, which will be handled through the corresponding pathways in the refined Bash pipeline.

### Sample Information

| Sample ID | SRA Accession | Tissue | Read Length | Layout |
|---|---|---|---|---|
| SKIN_1 | SRR7961333 | Skin | 50 bp | Single-end |
| SKIN_2 | SRR7961317 | Skin | 50 bp | Single-end |
| SKIN_3 | SRR7961304 | Skin | 50 bp | Single-end |
| SKIN_4 | SRR7961281 | Skin | 50 bp | Single-end |
| SKIN_5 | SRR33915754 | Skin | 75 bp | Single-end |
| COLON_1 | SRR7961292 | Colon | 50 bp | Single-end |
| COLON_2 | SRR7961262 | Colon | 50 bp | Single-end |
| COLON_3 | SRR7961278 | Colon | 50 bp | Single-end |
| COLON_4 | SRR7961254 | Colon | 50 bp | Single-end |
| COLON_5 | SRR2012208 | Colon | 76 bp | Paired-end |
| LUNG_1 | SRR7961222 | Lung | 50 bp | Single-end |
| LUNG_2 | SRR7961255 | Lung | 50 bp | Single-end |
| LUNG_3 | SRR7961291 | Lung | 50 bp | Single-end |
| LUNG_4 | SRR7961223 | Lung | 50 bp | Single-end |
| LUNG_5 | SRR7961221 | Lung | 50 bp | Single-end |

## Tools and Software

The refined 15-sample workflow is currently being developed locally in Ubuntu using the following bioinformatics and statistical software:

- **Falco v2.0.2** – quality assessment of raw, trimmed and host-unmapped RNA-seq reads.
- **Trimmomatic v0.41** – read preprocessing and quality trimming.
- **HISAT2 v2.2.3** – alignment to the human reference genome and extraction of host-unmapped reads.
- **Kraken2 v2.17.1** – taxonomic classification of host-unmapped reads.
- **Bracken v3.0.1** – species-level abundance estimation.
- **R v4.6.1** – downstream microbial community and statistical analysis.
- **vegan v2.7-5** – Bray-Curtis dissimilarity, PERMANOVA and PERMDISP analyses.
- **ggplot2 v4.0.3** – visualisation of microbial community patterns, including PCoA.



