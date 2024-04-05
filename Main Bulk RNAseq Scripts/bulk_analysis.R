library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(EnhancedVolcano)
library(factoextra)
library(rmarkdown)

set.seed(42)


############################################################################
# load data from count files and generate count matrix and sample metadata
############################################################################

counts=NULL
cfpattern = ".txt"
fs <- list.files(path=".", pattern=cfpattern)
cat("\n", length(fs), "count files detected")
cat("\nProcessing file:", fs[1])
counts <- read.table(file=fs[1], sep="\t", head=F, row.names=1)
colnames(counts) <- fs[1]

if (length(fs)>1) {
  for (i in 2:length(fs)) {
    cat("\nProcessing file:", fs[i])
    temp <- read.table(file=fs[i], sep="\t", head=F, row.names = 1)
    if(any(rownames(counts) != rownames(temp))) print(fs[i])
    counts[, fs[i]] <- temp[rownames(counts), 1]
  }
}
cat("\nCount matrix generated\n")

print(dim(counts))
print(head(counts))
colnames(counts) = gsub("\\.txt","",colnames(counts))
print(head(counts))

# generate sample metadata
samples = data.frame(sampleID=colnames(counts), 
                    condition=gsub("_[[:print:]]+$","",colnames(counts)))
samples$condition = factor(samples$condition, levels=c("NC", "BA"))
samples = samples[order(samples$condition), ]
rownames(samples) = samples$sampleID

counts = counts[, rownames(samples)]

saveRDS(counts, file="counts.rds")
saveRDS(samples, file="samples.rds")


############################################################################
# build DESeq2 model
############################################################################

PREFIX="BAvsNC"

colData = samples
colData$condition = factor(colData$condition, levels=c("NC","BA"))

counts.use = counts[, rownames(colData)]
counts.use = counts.use[which(rowSums(counts.use>1)>1), ] 

cat("Selected counts data:", dim(counts.use)[1], "genes x", dim(counts.use)[2], "samples")

dds <- DESeqDataSetFromMatrix(countData = counts.use, 
                               colData=colData, 
                               design = ~ condition)
dds <- DESeq(dds)

vsd <- vst(dds)
rld = rlog(dds)
print(assay(vsd)[1:3,1:3])
print(assay(rld)[1:3,1:3])

saveRDS(dds, file=paste0(PREFIX, ".dds.rds"))
saveRDS(vsd, file=paste0(PREFIX, ".vsd.rds"))
saveRDS(rld, file=paste0(PREFIX, ".rld.rds"))


############################################################################
# QC  
############################################################################

# PCA
g = plotPCA(vsd, intgroup=c("condition"), ntop=2000, returnData = F)
ggsave(filename=paste0(PREFIX, ".sample_pca.tiff"), plot = g, width=4.5, height=4, dpi=300, units="in", compression="lzw")


viz <- plotPCA(vsd, intgroup=c("condition"), ntop=2000, returnData = T)
g = ggplot(viz, aes(x=PC1, y=PC2, fill=condition, label=name)) + geom_point(size=4, shape=21, color="black") + scale_fill_manual(values=c("#006DDB","#FF6DB6"))
g = g + ggrepel::geom_text_repel(fill="white", ) + xlab("PC1: 52%") + ylab("PC2: 28%")
g = g + theme_bw() + theme(panel.grid = element_blank())
g
ggsave(filename=paste0(PREFIX, ".sample_pca.wlabel.png"), plot = g, width=4.5, height=3, dpi=300, units="in")


# Heatmap of sample correlations
tmp = t(scale(t(cor(assay(vsd)))))
colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
pheatmap(tmp, cluster_rows = T, cluster_cols = T, 
         #col=colors, 
         filename = paste0(PREFIX, ".sample_heatmap.png"), 
         width = 5, height = 4)

# Boxplot of gene expression distribution (after VST normalization) 
vsd.m <- melt(assay(vsd))
colnames(vsd.m) <- c("Gene","Sample","Expression")
g <- ggplot(vsd.m, aes(x=Sample, y=Expression)) 
g <- g + geom_boxplot(outlier.shape=NA)
g <- g + ylim(c(0, NA))
g <- g + theme_bw() + theme(panel.grid = element_blank())
g <- g + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
ggsave(filename = paste0(PREFIX, ".sample_boxplot.png"), width=4, height=3, dpi=300, units="in")


############################################################################
# get differential expression
############################################################################

res = results(dds)
summary(res)

ress = lfcShrink(dds, coef="condition_BA_vs_NC", res=res, type="apeglm")
summary(ress)

saveRDS(res, file=paste0(PREFIX, ".res.rds"))
saveRDS(ress, file=paste0(PREFIX, ".ress.rds"))

diff = subset(ress, abs(log2FoldChange) >= 1 & padj<0.05)
diff = diff[order(-diff$log2FoldChange), ]

write.table(diff, file=paste0(PREFIX, ".diff.shrinkFC2.q0.05.txt"), sep="\t", col.names=T, row.names=T, quote=FALSE)


############################################################################
# visualize differential expression 
############################################################################

# MA plot
plotMA(ress, xlab="Mean normalized RNA-seq signal")

# volcano plot
g = EnhancedVolcano(ress, title = NULL, subtitle = NULL,
                    lab = rownames(ress), 
                    selectLab = NULL,
                    FCcutoff = log2(2), 
                    pCutoff=0.05,
                    legendLabSize = 10,
                    legendPosition = "right",
                    x = 'log2FoldChange',
                    y = 'padj',
                    labSize = 1.2,
                    labFace = 'italic',
                    pointSize = 0.5,
                    axisLabSize = 9,
                    gridlines.major = FALSE,
                    gridlines.minor = FALSE,
                    caption = NULL,
                    borderWidth = 0.4)
ggsave(file=paste0(PREFIX, ".deg.volcano.tiff"), width=10, height=10, dpi=300, units="in", compression="lzw")

### heatmap of DEGs
mby.colors = colorRampPalette(c("#FF00FF", "#000000", "#FFFF00"))

col2 = colData
col2 = droplevels(col2)
col2$condition = factor(col2$condition, levels=c("NC","BA"))
ann_colors = list(condition=c(NC="#006DDB",BA="#FF6DB6"))

ann_col = data.frame(condition=col2[, "condition"])
rownames(ann_col) = rownames(col2)

diff$gene = rownames(diff)
mat = assay(vsd)[as.character(diff$gene), rownames(col2)]
mat = t(scale(t(mat), center = T, scale = T))
mat[mat>1.5] = 1.5
mat[mat<(-1.5)] = -1.5

pheatmap::pheatmap(as.matrix(mat), color = mby.colors(100), 
                   cluster_rows = F, cluster_cols = T, 
                   show_rownames = F,
                   show_colnames = T,
                   annotation_col = ann_col,
                   annotation_colors = ann_colors, border_color=NA,
                   filename = paste0(PREFIX, ".degs.shrinkfc", 2, "padj", 0.05, ".heatmap.tiff"), 
                   wdith=10, height=10)
