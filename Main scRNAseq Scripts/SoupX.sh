#BSUB -W 24:00
#BSUB -o SoupX.out
#BSUB -J SoupX
#BSUB -M 100000

module load R/4.1.1

Rscript SoupX.R