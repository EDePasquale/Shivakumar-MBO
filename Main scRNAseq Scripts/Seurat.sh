#BSUB -W 24:00
#BSUB -o Seurat.out
#BSUB -J Seurat
#BSUB -M 100000

module load R/4.1.1

Rscript Seurat.R