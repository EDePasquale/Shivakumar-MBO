#######################
#                     #
#  EMT and ECM Score  #
#     Correlation     #
#    Shivakumar-MBO   #
#    6 January 2025   #
#   Erica DePasquale  #
#                     #
#######################

# Load libraries
library(Seurat)
library(patchwork)
library(readxl)
library(ggplot2)
library(ggpubr)
library(plyr)

# Read in object
setwd("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5")
M<-readRDS("Seurat_8_MBO_withMeta_names_merged.rds")

# Bring in scoring function
source("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/200116_FunctionsGeneral.R")

# Make new folder for results
dir.create("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/ECM_EMT_Cor")
setwd("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/ECM_EMT_Cor")

# Read in ECM signauture first
signatures <- read_excel("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/ECM_signature/ECM_signature.xlsx", sheet = 1)
cexsize <- round(max(c(0.3, 0.5-nrow(M@assays[["RNA"]]@data)/40000)),2)
signatures <- lapply(as.list(signatures), function(y) y[!is.na(y)])
signatures <- lapply(signatures, intersect, rownames(M@assays[["RNA"]]@data))
CM.mean <- rowMeans(M@assays[["RNA"]]@data)
signScore_ECM <- lapply(names(signatures), function(g) {
  message(g)
  scoreSignature(CM = M@assays[["RNA"]]@data, signatures = signatures[[g]], CM.mean = CM.mean, verbose = T)
})
names(signScore_ECM) <- names(signatures)

# Then EMT signature
signatures <- read_excel("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/EMT_signatures/EMT_signatures.xlsx", sheet = 1)
cexsize <- round(max(c(0.3, 0.5-nrow(M@assays[["RNA"]]@data)/40000)),2)
signatures <- lapply(as.list(signatures), function(y) y[!is.na(y)])
signatures <- lapply(signatures, intersect, rownames(M@assays[["RNA"]]@data))
CM.mean <- rowMeans(M@assays[["RNA"]]@data)
signScore_EMT <- lapply(names(signatures), function(g) {
  message(g)
  scoreSignature(CM = M@assays[["RNA"]]@data, signatures = signatures[[g]], CM.mean = CM.mean, verbose = T)
})
names(signScore_EMT) <- names(signatures)

# Add scores to object
M@meta.data[["ECM"]]<-signScore_ECM$NABA_CORE_MATRISOME
M@meta.data[["EMT"]]<-signScore_EMT$EMTome_100

# Do correlation
cor.test(x=M$ECM, y=M$EMT)
# Pearson's product-moment correlation
# 
# data:  M$ECM and M$EMT
# t = 498.31, df = 27643, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9473985 0.9497605
# sample estimates:
#       cor 
# 0.9485927 

cor.test(x=M$ECM, y=M$EMT, method="spearman")
# Spearman's rank correlation rho
# 
# data:  M$ECM and M$EMT
# S = 2.3346e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.9336998 

# Result = very correlated!

# Make results tables for Hiro
M@meta.data[["merged_clusters_renamed"]]<-mapvalues(x=M$merged_clusters, from=unique(M$merged_clusters), to=c("Epithelial", "Mesenchymal", "Endothelial", "Intermediate 4", "Intermediate 3", "Intermediate 1", "Intermediate 2", "Intermediate 5", "Other"))
M$merged_clusters_renamed <- factor(M$merged_clusters_renamed, levels=c("Mesenchymal", "Epithelial", "Endothelial", "Intermediate 1",
                                                                        "Intermediate 2", "Intermediate 3", "Intermediate 4",
                                                                        "Intermediate 5", "Other"))

results_1=as.data.frame(cbind(Cell=Cells(M), ECM=M$ECM, EMT=M$EMT, Cluster=as.character(M$merged_clusters_renamed), Group=M$Diagnosis))
write.table(results_1, "data.txt", sep="\t", row.names=F, col.names=T, quote=F)

# Check for normal distribution
pdf(file = "distributions.pdf", width = 8.5, height = 7)
par(mar=c(2, 2, 2, 2))
  hist(M$ECM, col='steelblue', main='ECM')
  hist(M$EMT, col='steelblue', main='EMT')
dev.off()

# Do broken down stats for Hiro
Idents(M)<-"merged_clusters_renamed"
for(j in unique(M$Diagnosis)){
  print(j)
  for(i in unique(M$merged_clusters_renamed)){
    print(i)
    M_sub<-subset(M, idents=i)
    M_sub<-subset(M_sub, subset = Diagnosis == j)
    print(cor.test(x=M_sub$ECM, y=M_sub$EMT, method="spearman"))
  }
}

# Supplementary Figure 2D stats
M$Diagnosis <- factor(M$Diagnosis, levels=c("NC","BA"))
VlnPlot(M, features = "ECM", group.by = "merged_clusters_renamed", split.by = "Diagnosis", split.plot=T, alpha=0) + ggtitle("ECM Score in All Clusters") + ylab("ECM Score") + xlab(NULL)

Idents(M)<-"merged_clusters_renamed"

# BH (less strict)
results_ECM=as.data.frame(matrix(nrow=9, ncol=12))
iter=1
for(i in unique(M$merged_clusters_renamed)){
  print(i)
  M_sub<-subset(M, idents=i)
  results_2=as.data.frame(cbind(M_sub$ECM, as.character(M_sub$Diagnosis)))
  results_2[,1]<-as.numeric(results_2[,1])
  colnames(results_2)<-c("Value", "Group")
  X=results_2 %>% wilcox_test(Value ~ Group, detailed=T)
  results_ECM[iter, ]<-X
  iter=iter+1
}
colnames(results_ECM)=names(X)
row.names(results_ECM)=unique(M$merged_clusters_renamed)
results_ECM=cbind(results_ECM, p.adjust.BH=p.adjust(results_ECM$p, method="BH"))
write.table(results_ECM, "Wilcoxon_ECM_splitCluster_BH.txt")

results_EMT=as.data.frame(matrix(nrow=9, ncol=12))
iter=1
for(i in unique(M$merged_clusters_renamed)){
  print(i)
  M_sub<-subset(M, idents=i)
  results_2=as.data.frame(cbind(M_sub$EMT, as.character(M_sub$Diagnosis)))
  results_2[,1]<-as.numeric(results_2[,1])
  colnames(results_2)<-c("Value", "Group")
  X=results_2 %>% wilcox_test(Value ~ Group, detailed=T)
  results_EMT[iter, ]<-X
  iter=iter+1
}
colnames(results_EMT)=names(X)
row.names(results_EMT)=unique(M$merged_clusters_renamed)
results_EMT=cbind(results_EMT, p.adjust.BH=p.adjust(results_EMT$p, method="BH"))
write.table(results_EMT, "Wilcoxon_EMT_splitCluster_BH.txt")

# Bonferroni (more strict)
results_ECM=as.data.frame(matrix(nrow=9, ncol=12))
iter=1
for(i in unique(M$merged_clusters_renamed)){
  print(i)
  M_sub<-subset(M, idents=i)
  results_2=as.data.frame(cbind(M_sub$ECM, as.character(M_sub$Diagnosis)))
  results_2[,1]<-as.numeric(results_2[,1])
  colnames(results_2)<-c("Value", "Group")
  X=results_2 %>% wilcox_test(Value ~ Group, detailed=T)
  results_ECM[iter, ]<-X
  iter=iter+1
}
colnames(results_ECM)=names(X)
row.names(results_ECM)=unique(M$merged_clusters_renamed)
results_ECM=cbind(results_ECM, p.adjust.bonferroni=p.adjust(results_ECM$p, method="bonferroni"))
write.table(results_ECM, "Wilcoxon_ECM_splitCluster_bonferroni.txt")

results_EMT=as.data.frame(matrix(nrow=9, ncol=12))
iter=1
for(i in unique(M$merged_clusters_renamed)){
  print(i)
  M_sub<-subset(M, idents=i)
  results_2=as.data.frame(cbind(M_sub$EMT, as.character(M_sub$Diagnosis)))
  results_2[,1]<-as.numeric(results_2[,1])
  colnames(results_2)<-c("Value", "Group")
  X=results_2 %>% wilcox_test(Value ~ Group, detailed=T)
  results_EMT[iter, ]<-X
  iter=iter+1
}
colnames(results_EMT)=names(X)
row.names(results_EMT)=unique(M$merged_clusters_renamed)
results_EMT=cbind(results_EMT, p.adjust.bonferroni=p.adjust(results_EMT$p, method="bonferroni"))
write.table(results_EMT, "Wilcoxon_EMT_splitCluster_bonferroni.txt")


##############
# New requests May 2025

# "p-values (two-tailed Wilcoxon with BH correction) in ECM and EMT score

# Mesenchymal vs Endothelial
# Mesenchymal vs Epithelial
# Intermediate4 vs Endothelial
# Intermediate4 vs Mesenchymal
# Intermediate4 vs Epithelial
# 
# in NC vs NC, NC vs BA, BA vs BA"

M@meta.data[["merged_clusters_renamed"]]<-mapvalues(x=M$merged_clusters, from=unique(M$merged_clusters), to=c("Epithelial", "Mesenchymal", "Endothelial", "Intermediate 4", "Intermediate 3", "Intermediate 1", "Intermediate 2", "Intermediate 5", "Other"))
M$merged_clusters_renamed <- factor(M$merged_clusters_renamed, levels=c("Mesenchymal", "Epithelial", "Endothelial", "Intermediate 1",
                                                                        "Intermediate 2", "Intermediate 3", "Intermediate 4",
                                                                        "Intermediate 5", "Other"))
Idents(M)<-"merged_clusters_renamed"

test_conditions=as.data.frame(cbind(rbind("Mesenchymal_NC",
                                          "Mesenchymal_NC",
                                          "Intermediate 4_NC",
                                          "Intermediate 4_NC",
                                          "Intermediate 4_NC",
                                          "Mesenchymal_NC",
                                          "Mesenchymal_NC",
                                          "Intermediate 4_NC",
                                          "Intermediate 4_NC",
                                          "Intermediate 4_NC",
                                          "Mesenchymal_BA",
                                          "Mesenchymal_BA",
                                          "Intermediate 4_BA",
                                          "Intermediate 4_BA",
                                          "Intermediate 4_BA"),
                                    rbind("Endothelial_NC",
                                          "Epithelial_NC",
                                          "Endothelial_NC",
                                          "Mesenchymal_NC",
                                          "Epithelial_NC",
                                          "Endothelial_BA",
                                          "Epithelial_BA",
                                          "Endothelial_BA",
                                          "Mesenchymal_BA",
                                          "Epithelial_BA",
                                          "Endothelial_BA",
                                          "Epithelial_BA",
                                          "Endothelial_BA",
                                          "Mesenchymal_BA",
                                          "Epithelial_BA")))
colnames(test_conditions)<-c("Group1", "Group2")

myData=as.data.frame(cbind(M$ECM, as.character(M$Diagnosis), as.character(M$merged_clusters_renamed)))

results_ECM=as.data.frame(matrix(nrow=15, ncol=10))
results_ECM[,1:2]<-test_conditions
for(i in 1:nrow(test_conditions)){
  X=unlist(strsplit(test_conditions[i,1], split="_"))
  Y=unlist(strsplit(test_conditions[i,2], split="_"))
  myX=as.numeric(myData[intersect(which(myData$V3==X[1]),which(myData$V2==X[2])),1])
  myY=as.numeric(myData[intersect(which(myData$V3==Y[1]),which(myData$V2==Y[2])),1])
  if(length(myX>2) && length(myY>2)){
    results_ECM[i,3:8]<-unlist(wilcox.test(x=myX, y=myY, alternative="two.sided"))
  }else{
    print(paste0("Not enough observations: ", test_conditions[i,1], " ", length(myX>2), ", ", test_conditions[i,2], " ",length(myY>2)))
  }
  results_ECM[i,9]<-mean(myX)
  results_ECM[i,10]<-mean(myY)
}
colnames(results_ECM)=c("Group1", "Group2", names(unlist(wilcox.test(x=myX, y=myY, alternative="two.sided"))), "Mean Group1", "Mean Group2")
results_ECM=cbind(results_ECM, p.adjust.BH=p.adjust(results_ECM$p.value, method="BH"))
write.table(results_ECM, "Wilcoxon_ECM_May2025Requests_BH.txt")


myData=as.data.frame(cbind(M$EMT, as.character(M$Diagnosis), as.character(M$merged_clusters_renamed)))

results_EMT=as.data.frame(matrix(nrow=15, ncol=10))
results_EMT[,1:2]<-test_conditions
for(i in 1:nrow(test_conditions)){
  X=unlist(strsplit(test_conditions[i,1], split="_"))
  Y=unlist(strsplit(test_conditions[i,2], split="_"))
  myX=as.numeric(myData[intersect(which(myData$V3==X[1]),which(myData$V2==X[2])),1])
  myY=as.numeric(myData[intersect(which(myData$V3==Y[1]),which(myData$V2==Y[2])),1])
  if(length(myX>2) && length(myY>2)){
    results_EMT[i,3:8]<-unlist(wilcox.test(x=myX, y=myY, alternative="two.sided"))
  }else{
    print(paste0("Not enough observations: ", test_conditions[i,1], " ", length(myX>2), ", ", test_conditions[i,2], " ",length(myY>2)))
  }
  results_EMT[i,9]<-mean(myX)
  results_EMT[i,10]<-mean(myY)
}
colnames(results_EMT)=c("Group1", "Group2", names(unlist(wilcox.test(x=myX, y=myY, alternative="two.sided"))), "Mean Group1", "Mean Group2")
results_EMT=cbind(results_EMT, p.adjust.BH=p.adjust(results_EMT$p.value, method="BH"))
write.table(results_EMT, "Wilcoxon_EMT_May2025Requests_BH.txt")


