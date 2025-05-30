library(Seurat)
library(patchwork)
library(readxl)
library(ggplot2)
library(ggpubr)
library(plyr)

setwd("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5")
M<-readRDS("Seurat_8_MBO_withMeta_names_merged.rds")

source("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/200116_FunctionsGeneral.R")
message("colCustom()")
colCustom <- function(x, z=NULL, colors=c("#FFFFFF", "red")) {   # use grey to red as default
  if(is.null(z)) {   # just scale to min and max value
    #x <- scaleMinMax(x, min(x, na.rm=TRUE), max(x, na.rm=TRUE))
    x <- scaleMinMax(x, -0.1619832, 0.6405086) # CUSTOM FOR ECM SCORING BASED ON MAX and MAX VALUE of any group (so we can compare between conditions)
  } else if(length(z) == 1) {   # zscore
    x <- scaleMinMax(x, z = z, keepwithin = TRUE)
  } else {  # scale to min max as provided
    x <- scaleMinMax(x, min = z[1], max = z[2], keepwithin = TRUE)
  }

  m <- is.na(x)
  x[m] <- 0.5
  r <- colorRamp(colors)(x)
  y <- apply(r, 1, function(rr) rgb(rr[1], rr[2], rr[3], maxColorValue = 255))
  y[m] <- NA
  y
}

signatures <- read_excel("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/ECM_signature/ECM_signature.xlsx", sheet = 1)
cexsize <- round(max(c(0.3, 0.5-nrow(M@assays[["RNA"]]@data)/40000)),2)

dir.create("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/ECM_signature/Results")
setwd("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/ECM_signature/Results")

######
# NC #
######

M_NC<-subset(x = M, subset = Diagnosis == c("NC"))

# Color by signature score
message("\nCalculating signature scores")
# Signatures
signatures <- lapply(as.list(signatures), function(y) y[!is.na(y)])
signatures <- lapply(signatures, intersect, rownames(M_NC@assays[["RNA"]]@data))
# Average gene expression for scoreSignature
CM.mean <- rowMeans(M_NC@assays[["RNA"]]@data)
signScore <- lapply(names(signatures), function(g) {
  message(g)
  scoreSignature(CM = M_NC@assays[["RNA"]]@data, signatures = signatures[[g]], CM.mean = CM.mean, verbose = T)
})
names(signScore) <- names(signatures)

# Plot all signatures - UMAP
pdf(file = "UMAP_MBO_Seurat_ECM_scale_NC.pdf", width = 10, height = 10)
par(mar=c(4, 4, 4, 4))

for (n in names(signScore)) {
  mycol <- colItay(signScore[[n]])
  plot(M_NC@reductions[["umap"]]@cell.embeddings, pch = 16, cex = cexsize, col = mycol, main = n)
}
dev.off()

######
# BA #
######

M_BA<-subset(x = M, subset = Diagnosis == c("BA"))

# Color by signature score
message("Calculating signature scores")
# Signatures
signatures <- lapply(as.list(signatures), function(y) y[!is.na(y)])
signatures <- lapply(signatures, intersect, rownames(M_BA@assays[["RNA"]]@data))
# Average gene expression for scoreSignature
CM.mean <- rowMeans(M_BA@assays[["RNA"]]@data)
signScore <- lapply(names(signatures), function(g) {
  message(g)
  scoreSignature(CM = M_BA@assays[["RNA"]]@data, signatures = signatures[[g]], CM.mean = CM.mean, verbose = T)
})
names(signScore) <- names(signatures)

# Plot all signatures - UMAP
pdf(file = "UMAP_MBO_Seurat_ECM_scale_BA.pdf", width = 10, height = 10)
par(mar=c(4, 4, 4, 4))

for (n in names(signScore)) {
  mycol <- colItay(signScore[[n]])
  plot(M_BA@reductions[["umap"]]@cell.embeddings, pch = 16, cex = cexsize, col = mycol, main = n)
}
dev.off()

#######
# All #
#######

# Color by signature score
message("\nCalculating signature scores")
# Signatures
signatures <- lapply(as.list(signatures), function(y) y[!is.na(y)])
signatures <- lapply(signatures, intersect, rownames(M@assays[["RNA"]]@data))
# Average gene expression for scoreSignature
CM.mean <- rowMeans(M@assays[["RNA"]]@data)
signScore <- lapply(names(signatures), function(g) {
  message(g)
  scoreSignature(CM = M@assays[["RNA"]]@data, signatures = signatures[[g]], CM.mean = CM.mean, verbose = T)
})
names(signScore) <- names(signatures)
# Plot all signatures - UMAP
pdf(file = "UMAP_MBO_Seurat_ECM_scale.pdf", width = 10, height = 10)
par(mar=c(4, 4, 4, 4))

for (n in names(signScore)) {
  mycol <- colItay(signScore[[n]])
  plot(M@reductions[["umap"]]@cell.embeddings, pch = 16, cex = cexsize, col = mycol, main = n)
}
dev.off()

##########################
# Box Plot for ECM score #
##########################

M@meta.data[["ECM"]]<-signScore$NABA_CORE_MATRISOME
Idents(M)<-"Diagnosis"

NC=M@meta.data$ECM[which(M@meta.data$Diagnosis=="NC")]
BA=M@meta.data$ECM[which(M@meta.data$Diagnosis=="BA")]

myData=as.data.frame(rbind(cbind(NC,rep("NC",length(NC))),cbind(BA,rep("BA",length(BA)))))
colnames(myData)<-c("Score", "Group")
myData$Group=factor(myData$Group, levels = c("NC", "BA"))
myData$Score=as.numeric(myData$Score)

p <- ggplot(myData, aes(x=Group, y=Score, fill=Group)) + 
  geom_boxplot() +
  labs(title="ECM score",x="Group", y = "Score")

# uses Wilcoxon
pdf(file = "ECM_Score_Box_means.pdf", width = 5, height = 5)
par(mar=c(4, 4, 4, 4))
p + theme_classic() + stat_compare_means(label.y = 1.5)
dev.off()

# pdf(file = "ECM_Score_Violin_All_Clusters.pdf", width = 15, height = 5)
# par(mar=c(4, 4, 4, 4))
#   VlnPlot(M, features = "ECM", group.by = "cluster_names", split.by = "Diagnosis", split.plot=T, alpha=0) + ggtitle("ECM Score in All Clusters") + ylab("ECM Score") + xlab(NULL)
# dev.off()
# 
# pdf(file = "ECM_Score_Violin_Broad_Clusters.pdf", width = 5, height = 5)
# par(mar=c(4, 4, 4, 4))
#   VlnPlot(M, features = "ECM", group.by = "cluster_names_redu", split.by = "Diagnosis", split.plot=T, alpha=0) + ggtitle("ECM Score in Broad Clusters") + ylab("ECM Score") + xlab(NULL)
# dev.off()

M@meta.data[["merged_clusters_renamed"]]<-mapvalues(x=M$merged_clusters, from=unique(M$merged_clusters), to=c("Epithelial", "Mesenchymal", "Endothelial", "Intermediate 4", "Intermediate 3", "Intermediate 1", "Intermediate 2", "Intermediate 5", "Other"))
M$merged_clusters_renamed <- factor(M$merged_clusters_renamed, levels=c("Mesenchymal", "Epithelial", "Endothelial", "Intermediate 1",
                                                        "Intermediate 2", "Intermediate 3", "Intermediate 4",
                                                        "Intermediate 5", "Other"))

M$Diagnosis <- factor(M$Diagnosis, levels=c("NC","BA"))
pdf(file = "ECM_Score_Violin_All_Clusters.pdf", width = 15, height = 5)
par(mar=c(4, 4, 4, 4))
VlnPlot(M, features = "ECM", group.by = "merged_clusters_renamed", split.by = "Diagnosis", split.plot=T, alpha=0) + ggtitle("ECM Score in All Clusters") + ylab("ECM Score") + xlab(NULL)
dev.off()


########
#install.packages("binovisualfields")
library(binovisualfields)

colors=c("#660220", "#b01b2f", "#d46151", "#f2a585", "#fcdbc8", "#f7f7f7", "#d2e5ef", "#94c5dd", "#4794c1", "#2668aa", "#083160")

pdf(file = "ECM_Legend.pdf", width = 10, height = 5)
par(mar=c(4, 4, 4, 4))
  layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
  plot(1:11, 1:11, pch = 19, cex=2, col = colors)

  legend_image <- as.raster(matrix(colors, ncol=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'ECM Score')
  #text(x=1.5, y = seq(0,1,l=5), labels = seq(round(min(signScore[[1]]),3),round(max(signScore[[1]]),3),l=5)) 
  #text(x=1.5, y = seq(0,1,l=5), labels = seq(round(-0.1619832,3),round(0.6405086,3),l=5)) # <- set true min and mix (was set by all only before, when NA ranges were slightly larger)
  text(x=1.5, y = seq(0,1,l=5), labels = c(-0.162, 0.039, 0.24, 0.440, 0.641)) # <- got weird visual effects with this so I did it manually
  rasterImage(legend_image, 0, 0, 1,1)
dev.off()

