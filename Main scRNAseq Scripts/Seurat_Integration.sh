#BSUB -W 24:00
#BSUB -o Seurat_Integration_0.5.out
#BSUB -J Seurat_Integration_0.5
#BSUB -M 100000

module load R/4.1.1

Rscript Seurat_Integration.R