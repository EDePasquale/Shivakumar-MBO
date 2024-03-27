##############################
#                            #
#     Prepare files for      #
#        CellphoneDB         #
#      Shivakumar - MBO      #
#     Erica DePasquale       #
#        11 May 2023         #
#                            #
##############################

library(Seurat)
library(SeuratObject)
library(Matrix)

# From: https://github.com/ventolab/CellphoneDB/blob/master/notebooks/0_prepare_your_data_from_Seurat.ipynb

# 0. Make working directory
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5/")
dir.create("CellPhoneDB_BA")
setwd("CellPhoneDB_BA")
dir.create("Shivakumar_MBO_counts_mtx")

# 1. Load seurat object
so <- readRDS("/Volumes/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5/Seurat_8_MBO_withMeta_names_merged.rds")
Idents(so) = so$merged_clusters
so<-subset(so, subset = Diagnosis == "BA")

# 2. Write gene expression in mtx format
writeMM(so@assays$RNA@counts, file = 'Shivakumar_MBO_counts_mtx/matrix.mtx')
write(x = rownames(so@assays$RNA@counts), file = "Shivakumar_MBO_counts_mtx/features.tsv")
write(x = colnames(so@assays$RNA@counts), file = "Shivakumar_MBO_counts_mtx/barcodes.tsv")

# 3. Generate your meta
#In this example, our input is an anndata containing the cluster/celltype information in metadat$'cell_type'
#The object also has metadat$'lineage' information wich will be used below for a hierarchical DEGs approach.
so@meta.data$Cell = rownames(so@meta.data)
df = so@meta.data[, c('Cell', 'merged_clusters')]
write.table(df, file ='Shivakumar_MBO_meta.tsv', sep = '\t', quote = F, row.names = F)

# 4. Compute DEGs (optional)
## OPTION 1 - compute DEGs for all cell types
DEGs <- FindAllMarkers(so,
                       test.use = 'LR',
                       verbose = T,
                       only.pos = T,
                       random.seed = 1,
                       logfc.threshold = 0.2,
                       min.pct = 0.1,
                       return.thresh = 0.05)

fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_log2FC > 0.1)
fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')] 
write.table(fDEGs, file ='Shivakumar_MBO_DEGs.tsv', sep = '\t', quote = F, row.names = F)

# 5. Run cellphoneDB (not in R)
# cellphonedb method degs_analysis  \
#   Shivakumar_MBO_meta.tsv  \
#   Shivakumar_MBO_counts_mtx  \
#   Shivakumar_MBO_DEGs.tsv  \
#   --microenvs Shivakumar_MBO_microenviroments.tsv  \ #optional
#   --counts-data hgnc_symbol  \
#   --database database/database/cellphonedb_user_2021-06-29-11_41.db \
#   --threshold 0.1


