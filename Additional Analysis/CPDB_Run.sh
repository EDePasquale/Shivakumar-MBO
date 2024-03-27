#BSUB -W 24:00
#BSUB -o Shivakumar_BA_full.out
#BSUB -J Shivakumar_BA_full
#BSUB -M 100000

module load python3
source ~/cpdb/bin/activate
module load R/4.1.1

cellphonedb method degs_analysis  \
    /data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5/CellPhoneDB_BA/Shivakumar_MBO_meta.tsv  \
    /data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5/CellPhoneDB_BA/Shivakumar_MBO_counts_mtx  \
    /data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5/CellPhoneDB_BA/Shivakumar_MBO_DEGs.tsv  \
    --output-path /data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5/CellPhoneDB_BA/Shivakumar_DEG_results\
    --verbose \
    --counts-data hgnc_symbol

#bsub < /data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Scripts_used/Shivakumar_MBO/CPDB_Run.sh