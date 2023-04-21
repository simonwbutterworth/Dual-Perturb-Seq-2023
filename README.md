# Dual-Perturb-Seq-2023

Code repository to reproduce single-cell perturb-seq data analysis in the paper:

>**Mapping host-microbe transcriptional interactions by dual perturb-seq**
<br>Simon Butterworth, Kristina Kordova, Sambamurthy Chandrasekaran, Kaitlin K Thomas, Francesca Torelli, Eloise J Lockyer, Anita Koshy, Moritz Treeck<br>

## Introduction

This set of experiments and analyses use CRISPR-Cas9 knockout screening in intracellular *Toxoplasma gondii* parasites combined with single-cell RNA-sequencing of infected host cells to identify parasite effector proteins that modify transcription in the host cell, a method we term dual perturb-seq. Three dual perturb-seq experiments are described in the paper:
- A pilot experiment using *T. gondii* parasites transfected with two perturb-seq plasmids, targeting the UPRT and MYR1 genes.
- A larger validation experiment using parasites transfected with a panel of 24 perturb-seq plasmids.
- A screen of 256 *T. gondii* genes encoding putative secreted proteins.

*T. gondii* CRISPR mutant-infected human foreskin fibroblasts were sequenced using the 10x Genomics single-cell 3' gene expresion kit with Feature Barcoding for CRISPR screening. For each of the first two experiments, a single Chromium Controller channel was loaded. For the screen, three replicates were carried out in unstimulated HFFs and two in interferon-gamma-stimulated HFFs, with 2-4 channels per replicate, for a total of 16 channels. Demultiplexed FASTQ files have been deposited at the Gene Expression Omnibus under the accession number **GSE229505**. 

Pre-processing of FASTQ files to generate count matrices was carried out using Cell Ranger v6.1.2 (10x Genomics) with a custom multi-species reference genome comprising *Homo sapiens* GRCh38 Ensembl release 95 and *Toxoplasma gondii ME49* TGA4 Ensembl Protists release 42, plus a custom reference of sgRNA protospacer sequences created with Cell Ranger `mkref`. This pre-processing is not included in this repository, yet can be straightforwardly replicated with FASTQ files from GEO and the `cellranger_sgrna_reference` files included here, using Cell Ranger `count` with the --libraries flag to specify samples as "Gene Expression" or "CRISPR Guide Capture," the --feature-ref flag to pass the sgRNA protospacer reference files, and the --include-introns flag to retain intronic reads. All other settings are left as default.

## Analysis overview

Data analysis in this repository starts from the Cell Ranger "filtered_feature_bc_matrix" outputs, which should be downloaded from GEO into a `data` directory prior to running these scripts. All analyses are carried out in R using Seurat v4.2.0.

Refer to Figure S1 of Butterworth *et al.* 2023 for an overview of the data analysis strategy described in this repository.

## Repository contents

- `inputs`
  - `cellranger_sgrna_reference` - Cell Ranger --feature-ref inputs (used for pre-processing with Cell Ranger)
  - `seurat_sgrna_reference` - Lookup tables to assign sgRNA target identities to cell barcodes in Seurat object metadata
  - `gene_sets` - Pathway Interaction Database gene sets and custom gene sets for *T. gondii* effector-regulated genes in .gmt format
  - `phenotype_scores` - *In vitro* CRISPR phenotype scores from Sidik *et al.* 2016 for genes included in the dual perturb-seq screen 
- `script` - .R files with analysis code. FUNCTIONS.R contains libraries and functions used by all other scripts
- `outputs` - Outputs from script files
