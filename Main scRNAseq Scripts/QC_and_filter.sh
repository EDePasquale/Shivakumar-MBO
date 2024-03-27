#BSUB -W 24:00
#BSUB -o QC_and_filter.out
#BSUB -J QC_and_filter
#BSUB -M 100000

module load R/4.1.1

Rscript QC_and_filter.R