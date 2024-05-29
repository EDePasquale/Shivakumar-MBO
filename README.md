# Shivakumar-MBO

## Paper Citation
[to be determined]

## Stable Code Release
Release 1.0.0: [![DOI](https://zenodo.org/badge/777389401.svg)](https://zenodo.org/doi/10.5281/zenodo.10932740)

Release 1.0.1: [![DOI](https://zenodo.org/badge/777389401.svg)](https://zenodo.org/doi/10.5281/zenodo.11387118)

## Project Description
In this project, the authors investigate the pathogenesis of biliary atresia (BA) by creating and sequencing patient-derived multilineage biliary organoids (MBO). Bulk RNA sequencing and single-cell RNA sequencing were used to investigate markers of epithelial-mesenchymal transition (EMT) in these MBOs. 

For the bulk RNAseq, we performed quality assessment and pre-processing using FASTQC, Trim Galore, Cutadapt, and SAMtools. We aligned the reads to hg38 genome using Bowtie2 and removed low-quality alignments and PCR duplicates using SAMtools and Picard MarkDuplicates. We quantified gene expression using htseq-count, analyzed differential expression using DESeq2, and performed PCA using FactoMineR. PCAs, box plots, heatmaps and volcano plots were visualized using ggplot2, pheatmap and EnhancedVolcano.

For the proteomics data, we analyzed differential expression using stats and visualized heatmap using pheatmap.

For the single-cell portion of the paper, we used Cell Ranger to generate gene-by-cell expression matrices, SoupX to correct for ambient RNAs, and Seurat to perform QC and data integration, as well as most of the plots and statistics. In addition to generating the integrated dataset, we quantified differential cell type proportions between BA MBOs and normal control (NC) MBOs with Seurat, performed gene set enrichment using Enrichr, scored cells for known EMT transcriptional markers, and used CellPhoneDB to predict cell-cell communication.

## Processing Steps
### Bulk RNA Sequencing Scripts

#### Main Bulk RNAseq Pipeline
The following script was used to process all bulk RNAseq data and to generate figures:
1. All analyses and figure generation (bulk_analysis.R)

### Proteomics Scripts
The following files include codes for specific analyses on proteomics:
1. Proteomics_analysis_code.R
2. heatmap_code.R

### Single Cell RNA Sequencing Scripts

#### Main scRNAseq Pipeline
The following scripts were used to perform the following basic tasks to generate the integrated dataset:
1. Cell Ranger (cellranger6.sh, run.sh)
2. SoupX (SoupX.sh and .R)
3. QC and filtering (QC_and_Filtering.sh and .R)
4. Make Seurat objects (Seurat.sh and .R)
5. Sample integration (Seurat_Integration.sh and .R)

#### Additional Analysis
The following scripts were used to perform specific analyses on the integrated dataset:
1. Differential cell type proportions (UMAP_and_Differential_Proportions.R)
2. Gene set enrichment analysis (Enrichr.R)
3. EMT scoring and statistics (EMT_Scoring.R)
4. Cell-cell communication prediction with CellPhoneDB (CPDB_Prep.R, CPDB_Run.sh, CPDB_Post.R)

#### Accessory Files and Scripts
The following scripts and files are needed for running some of the above scripts. These include:
1. Gene information metadata sheet, needed for QC and filtering (gene_info_GRCh38_PvG200204.txt)
2. Signature genes for EMT scoring. The set described in the paper is the second column (EMT_signatures.xlsx)
3. Functions used in gene scoring, prior to modification, from Single-cell_BPDCN_Functions.R at https://github.com/petervangalen/Single-cell_BPDCN/ (200116_FunctionsGeneral.R)

##### Note: All additional figure panels were generated outside R (e.g., in Prism) and the code is not included in this repository.
