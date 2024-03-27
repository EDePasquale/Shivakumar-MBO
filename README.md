# Shivakumar-MBO
Code for analysis and figure generation for Ayabe et al. 2024

## Paper Citation
[to be determined]

## Project Description
In this project, the authors 

## Processing Steps
### Bulk RNA Sequencing Scripts

### Single Cell RNA Sequencing Scripts

#### Main scRNAseq Pipeline
The following scripts were used to perform the following basic tasks to generate the integrated dataset:
1. CellRanger (cellranger6.sh, run.sh)
2. SoupX (SoupX.sh and .R)
3. QC and filtering (QC_and_Filtering.sh and .R)
4. Make Seurat objects (Seurat.sh and .R)
5. Sample integration (Seurat_Integration.sh and .R)

#### Additional Analysis
The following scripts were used to perform specific analyses on the integrated dataset:
1. Differential cell type proportions ()
2. Gene set enrichment analysis ()
3. EMT scoring ()
4. Cell-cell communication prediction with CellPhoneDB ()

#### Accessory Files and Scripts
The following scripts and files are needed for running some of the above scripts. These include:
1. Gene information metadata sheet, needed for QC and filtering (gene_info_GRCh38_PvG200204.txt)

