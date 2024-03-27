library(Seurat)
#devtools::install_github('satijalab/seurat-data')
#library(SeuratData)
library(patchwork)
library(gdata)

setwd("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO")

# Load objects
CO15_MBO_NC=readRDS("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/QC/CO15-MBO-NC-20210816/CO15-MBO-NC-20210816_filtered_matrices/Seurat.rds")
CO36_MBO_NC=readRDS("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/QC/CO36-MBO-NC-20211103/CO36-MBO-NC-20211103_filtered_matrices/Seurat.rds")
CO53_MBO_NC=readRDS("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/QC/CO53-MBO-NC-20211201/CO53-MBO-NC-20211201_filtered_matrices/Seurat.rds")
CO55_MBO_NC=readRDS("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/QC/CO55-MBO-NC-20211027/CO55-MBO-NC-20211027_filtered_matrices/Seurat.rds")
CO21_MBO_BA=readRDS("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/QC/CO21_MBO-BA-20220621/CO21_MBO-BA-20220621_filtered_matrices/Seurat.rds")
CO22_MBO_BA=readRDS("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/QC/CO22-MBO-BA-20211117/CO22-MBO-BA-20211117_filtered_matrices/Seurat.rds")
CO68_MBO_BA=readRDS("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/QC/CO68-MBO-BA-20220302/CO68-MBO-BA-20220302_filtered_matrices/Seurat.rds")
CO72_MBO_BA=readRDS("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/QC/CO72-MBO-BA-20220509/CO72-MBO-BA-20220509_filtered_matrices/Seurat.rds")


ifnb.list <- list(CO15_MBO_NC, CO36_MBO_NC, CO53_MBO_NC, CO55_MBO_NC, CO21_MBO_BA, CO22_MBO_BA, CO68_MBO_BA, CO72_MBO_BA)
remove(CO15_MBO_NC, CO36_MBO_NC, CO53_MBO_NC, CO55_MBO_NC, CO21_MBO_BA, CO22_MBO_BA, CO68_MBO_BA, CO72_MBO_BA)

resolution=0.5
dir.create(paste0("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_", resolution))

for(i in 1:length(ifnb.list)){
  ifnb.list[[i]] <- NormalizeData(ifnb.list[[i]])
  ifnb.list[[i]] <- FindVariableFeatures(ifnb.list[[i]], selection.method = "vst", nfeatures = 2000)
  #ifnb.list[[i]]@assays[["RNA"]]@var.features=setdiff(ifnb.list[[i]]@assays[["RNA"]]@var.features, cc.genes) # for cell cycle removal in anchors
}

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# this command creates an 'integrated' data assay
immune.combined2 <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined2) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined2 <- ScaleData(immune.combined2, verbose = FALSE)
immune.combined2 <- RunPCA(immune.combined2, npcs = 30, verbose = FALSE)
immune.combined2 <- RunUMAP(immune.combined2, reduction = "pca", dims = 1:30)
immune.combined2 <- RunTSNE(immune.combined2, reduction = "pca", dims = 1:30)
immune.combined2 <- FindNeighbors(immune.combined2, reduction = "pca", dims = 1:30)
immune.combined2 <- FindClusters(immune.combined2, resolution = as.numeric(resolution))

# Visualization
p1 <- DimPlot(immune.combined2, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(immune.combined2, reduction = "umap", label = TRUE, repel = TRUE)

pdf(file = paste0("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_", resolution, "/UMAP_sample_and_cluster.pdf"), width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p1 + p2
dev.off()

pdf(file = paste0("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_", resolution, "/UMAP_split_by_sample.pdf"), width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(immune.combined2, reduction = "umap", split.by = "orig.ident")
dev.off()

# Visualization
p1 <- DimPlot(immune.combined2, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(immune.combined2, reduction = "tsne", label = TRUE, repel = TRUE)

pdf(file = paste0("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_", resolution, "/TSNE_sample_and_cluster.pdf"), width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p1 + p2
dev.off()

pdf(file = paste0("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_", resolution, "/TSNE_split_by_sample.pdf"), width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(immune.combined2, reduction = "tsne", split.by = "orig.ident")
dev.off()

#############

source("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/200116_FunctionsGeneral.R")

signatures <- read.xls("/data/GI-Informatics/DePasquale/Signature_Files/220519_signatures_top150.xlsx", header = T, sheet = 1)
cexsize <- round(max(c(0.3, 0.5-nrow(immune.combined2@assays[["RNA"]]@data)/40000)),2)

# Color by signature score
message("\nCalculating signature scores")
# Signatures
signatures <- lapply(as.list(signatures), function(y) y[!is.na(y)])
signatures <- lapply(signatures, intersect, rownames(immune.combined2@assays[["RNA"]]@data))
# Average gene expression for scoreSignature
CM.mean <- rowMeans(as.matrix(immune.combined2@assays[["RNA"]]@data))
signScore <- lapply(names(signatures), function(g) {
  message(g)
  scoreSignature(CM = as.matrix(immune.combined2@assays[["RNA"]]@data), signatures = signatures[[g]], CM.mean = CM.mean, verbose = T)
})
names(signScore) <- names(signatures)

# Plot all signatures
pdf(file = paste0("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_", resolution, "/UMAP_8_MBO_Seurat.pdf"), width = 10, height = 10)
par(mar=c(4, 4, 4, 4))

for (n in names(signScore)) {
  mycol <- colItay(signScore[[n]])
  plot(immune.combined2@reductions[["umap"]]@cell.embeddings, pch = 16, cex = cexsize, col = mycol, main = n)
}

dev.off()

# Plot all signatures
pdf(file = paste0("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_", resolution, "/TSNE_8_MBO_Seurat.pdf"), width = 10, height = 10)
par(mar=c(4, 4, 4, 4))

for (n in names(signScore)) {
  mycol <- colItay(signScore[[n]])
  plot(immune.combined2@reductions[["tsne"]]@cell.embeddings, pch = 16, cex = cexsize, col = mycol, main = n)
}

dev.off()

# Plot Seurat UMAP
pdf(file = paste0("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_", resolution, "/UMAP_8_MBO_Seurat_Groups.pdf"), width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = immune.combined2, reduction = "umap",label= TRUE)
dev.off()

# Plot Seurat TSNE
pdf(file = paste0("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_", resolution, "/TSNE_8_MBO_Seurat_Groups.pdf"), width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = immune.combined2, reduction = "tsne",label= TRUE)
dev.off()


saveRDS(immune.combined2, file = paste0("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_", resolution, "/Seurat_8_MBO.rds"))

temp=as.data.frame(cbind(immune.combined2@meta.data[["orig.ident"]],immune.combined2@meta.data[[paste0("integrated_snn_res.", resolution)]]))
write.table(table(temp), paste0("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_", resolution, "/Seurat_8_MBO_table.txt"), sep="\t", quote=F)
























# 
# 
# #####
# immune.combined=readRDS("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem_tcell/T_cell_0.5/Seurat.rds")
# immune.combined2=readRDS("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/Seurat_int_0.5.rds")
# 
# full_idents=as.data.frame(immune.combined2@active.ident)
# temp=cbind(immune.combined@active.ident, as.character(full_idents[names(immune.combined@active.ident),]))
# clus_names=as.data.frame(cbind(0:15, c("CD8T", "B", "MidEryth", "CTL", "Mono", "NK","LateEryth", "EarlyEryth", "CD4T", "PreB",
#                                        "ProB", "HSCProg", "cDC", "ncMono", "Plasma",  "pDC")))
# temp2=as.data.frame(temp)
# temp3=left_join(temp2, clus_names, by=c("V2"="V1"))
# immune.combined@meta.data[["old.ident"]]=as.factor(temp3[,3])
# 
# # Plot Seurat UMAP
# pdf(file = paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem_tcell/T_cell_0.5/UMAP_Full_Seurat_Groups.pdf"), width = 6, height = 6)
# par(mar=c(2, 2, 2, 2))
# DimPlot(object = immune.combined, reduction = "umap",label= TRUE, group.by = "old.ident")
# dev.off()
# #saveRDS(immune.combined, file = paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_", resolution, "_FLT3/Seurat_int_", resolution,".rds"))
# 
# 
# 
# 
# 
# 
# 
# #####
# #overlay ICGS-merge labels on Seurat generated UMAP plot for comparison
# immune.combined2=readRDS("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5/Seurat_11_Bile.rds")
# immune.combined2=readRDS("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem_tcell_2/Seurat_int_0.5.rds")
# immune.combined2=readRDS("/Users/erica/Documents/Projects/Data_files_temp/Donor/BPDCN712/Seurat.rds")
# #immune.combined2=readRDS("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem_tcell/T_cell_0.1/Seurat.rds")
# #ICGS=read.table("~/Documents/Projects/Data_files_temp/Donor/ICGS_Merge/ICGS-NMF/FinalGroups.txt", sep="\t")
# #ICGS=cbind(ICGS[,1], ICGS)
# #ICGS[,1]=gsub(".*[.]", "", ICGS[,1])
# #ICGS[,2]=gsub("[.].*", "", ICGS[,2])
# Seu=immune.combined2@reductions[["umap"]]@cell.embeddings
# Seu=cbind(row.names(Seu), Seu)
# library(dplyr)
# Seu=cbind(Seu, Seu[,1])
#Seu[,1]=gsub("_[0-9]$", "", Seu[,1])
# #Seu_ICGS=left_join(as.data.frame(Seu), as.data.frame(ICGS))
# Seu_ICGS=as.data.frame(Seu)
# 
# #start plotting
# library(randomcoloR)
# Seu_ICGS=Seu_ICGS[complete.cases(Seu_ICGS), ] #remove rows with NA
# n <- length(unique(Seu_ICGS$V4)) #how many colors do we need?
# palette <- distinctColorPalette(n) #make them distinct
# pie(rep(1,n), col=palette) #make a pie chart so you can confirm they are distinct
# palette_clus=as.data.frame(cbind(unique(Seu_ICGS$V4), palette)) #left_join doesn't work on row names, so make them a column
# 
# palette_clus[,1]=as.numeric(palette_clus[,1]) #make it numeric because factors suck
# 
# Seu_ICGS=left_join(Seu_ICGS, palette_clus, by=c("V2"="V1")) #attach these to the data.frame so we can get colors per cell
# 
# 
# Seu_ICGS_plot = as.matrix(Seu_ICGS[, 2:3]) #pull out the columns of coords we need to plot
# 
# pdf(file = paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_", resolution, "/UMAP_ICGS_Groups.pdf"), width = 6, height = 6)
# plot(Seu_ICGS_plot, pch = 16, cex=0.4, col = as.character(Seu_ICGS$palette), main = "ICGS clusters on Seurat UMAP") #plot with the wrong colors???
# legend("topleft", legend=unique(Seu_ICGS$V2), cex=0.5, pch=16, col=as.character(unique(Seu_ICGS$palette))) #legend
# dev.off()
# 
# #color plot by gene
# for(i in gene_list){
#   gene=i
#   # SELL, CCR7, PTPRC
#   #filt_exp=immune.combined2@assays[["RNA"]]@data[,Seu_ICGS[,4]]
#   #filt_exp=temp2
#   filt_exp=immune.combined2@assays[["RNA"]]@data
#   ind=which(row.names(filt_exp)==gene)
#   if(length(ind)!=0){
#     #df=data.frame(cell=Seu_ICGS[,1], exp=filt_exp[ind,])
#     df=data.frame(cell=Seu_ICGS[,1], exp=filt_exp[ind,])
#     rbPal <- colorRampPalette(c('grey','red'))
#     df$Col <- rbPal(10)[as.numeric(cut(df$exp,breaks = 10))]
#     pdf(paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_", resolution, "_FLT3/", gene,  ".pdf"), width=6, height=6)
#     plot(Seu_ICGS_plot, pch = 16, cex=0.2, col = df$Col, main = gene)
#     dev.off()
#   }else{
#     print(paste0("Gene ", gene, " not found"))
#   }
# }
# 
# 
# #cant find: CD303 (BDCA2), CD304 (BDCA4), CD123 (IL-3R), CD45RA
# 
#make PDF plot
#resolution=0.3
#gene_list=c("CD34", "GATA2", "MEIS1", "FCER1A", "CD14", "LYZ", "FCGR3A", "MS4A7", "MS4A1", "GNLY", "NKG7", "NCAM1", "CD8A", "IL7R", "CCR7", "KLF1")
#gene_list=c("IL7R", "PDCD1", "CXCR5", "IFNG", "PTPRC", "CCR7", "SELL", "CD4", "CD8A")
#gene_list=c("TNF", "IL2", "TBX21", "CCR4", "IL10", "IRF4", "CCR6", "KLRB1", "RORC")
# gene_list=c("PTPRC", "CD3", "CD4", "CD8A", "CD19", "MS4A1", "FCGR3A", "NCAM1", "IL3RA",
#             "HLA-DRA", "IL7R", "CD14", "CD2", "KLRK1", "MUC1", "KLRC1", "B3GAT1", "CD27",
#             "CCR7", "CD28", "FAS", "PDCD1", "CD38", "IL2RA", "ENTPD1", "CD24", "CCR5",
#             "CCR6", "CXCR3", "CXCR5", "ITGAX", "THBD", "CD1C")
#gene_list=c("CD8A", "CD4", "ITGB1", "CCR7", "CD27", "IL7R")
# # gene_list=c("CD68", "CD3D", "FCGR3A", "MS4A1", "SDC1", "PECAM1", "ANXA4", "SOX9", "KRT7", "S100A12", "S100A8", "CD4", "CD8A")
# 
# #Wang NK markers
# # gene_list=c("GNLY", "FGFBP2", "FCGR3A", "GZMB", "TYROBP", #GNLY_NK
#             "PTGDS", "FCER1G", "MYOM2", "AREG", "SPON2", #PTGDS_NK
#             "XCL2", "CMC1", "KLRC3", "KLRF1", "TYROBP", #XCL2_NK
#             "XCL1", "FGER1G", "XCL2", "AREG", "IL2RB") #XCL1_NK
# 
# #Wang T markers
# gene_list=c("CD8B", "GZMK", "RGS1", "CD8A", "COTL1", #CD8B_CD8
#             "CCL20", "KLRB1", "TRAV1-2", "1L7R", "SLC4A10", #CCL20_CD8
#             "IL7R", "GPR183", "LTB", "RGCC", "LEF1", #IL7R_CD4
#             "TRBV9", "TRAV38-2DV8", "TRBV13", "RP11-291B21.2", "TRGV5", #TRBV9_CD8
#             "TRDV2", "TRGV9", "KLRC1", "TRDC", "TRAC") #TRDV2_ydT
# 
# #Wang B Markers 
# gene_list=c("TCL1A", "TXNIP", "BTG1", "CD37", "FCER2",
#             "GPR183", "TNF", "COTL1", "AC079767", "NR4A2",
#             "MT2A", "FGR", "FCRL3", "FCRL5", "CIB1")
# 
# #Wang plasma markers
# gene_list=c("IGKVID-16", "IGKV3D-20", "IGKV1-16", "IGHV3-72", "IGKV2D-28",
#             "HIST1H4C", "RRM2", "GAPDH", "HMGB2", "IGHV6-1",
#             "IGLV3-1", "IGKV1D-39", "CH17-224D4.2", "IGKV1D-12", "IGLV3-21")
# 
# #gene_list=c("ITGB1", "ITGA4", "ITGA5", "ITGAL", "ITGAM", "ITGB2", "CD2", "CD44", "ICAM1", "CD58")
# #gene_list=c("HAVCR2", "PDCD1", "LAG3", "ENTPD1", "CTLA4", "TIGIT", "TNFRSF9", "CD27", "MYO7A", "WARS", "CXCL13", "TOX", "LAYN", "PHLDA1", "SNAP47", "TCF1", "CXCR5", "FOXP3")
# #pdf(file = paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_", resolution, "_FLT3/UMAP_TCELL_Markers2.pdf"), width = 6, height = 6)
# #pdf(file = paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_",resolution,"_ccrem_YandO/UMAP_basic_Tcell_Markers.pdf"), width = 5, height = 5)
# pdf(file = paste0("/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5/UMAP_Markers.pdf"), width = 7.5, height = 5)
# #pdf(file = paste0("~/Documents/Projects/Data_files_temp/Donor/BPDCN712/UMAP_exhaust_Tcell_Markers.pdf"), width = 10, height = 5)
# par(mar=c(1, 1, 1, 1))
# layout(matrix(1:6, ncol=3, nrow=2, byrow=TRUE))
# 
# for(n in 1:length(gene_list)){
#   gene=gene_list[n]
#   #filt_exp=immune.combined2@assays[["RNA"]]@data[,Seu_ICGS[,4]]
#   filt_exp=immune.combined2@assays[["RNA"]]@data
#   ind=which(row.names(filt_exp)==gene)
#   if(length(ind)!=0){
#     df=data.frame(cell=Seu_ICGS[,1], exp=filt_exp[ind,])
#     rbPal <- colorRampPalette(c('grey','red'))
#     df$Col <- rbPal(10)[as.numeric(cut(df$exp,breaks = 10))]
#     plot(Seu_ICGS[,2:3], pch = 16, cex=0.2, col = df$Col, main = gene, axes=FALSE, frame.plot=TRUE)
#     Axis(side=1, labels=FALSE, tick=F)
#     Axis(side=2, labels=FALSE, tick=F)
#   }else{
#     print(paste0("Gene ", gene, " not found"))
#   }
# }
# 
# dev.off()
# 
# 
# 
# 
# #ignore
# # setwd("~/Documents/Projects/Data_files_temp/Donor/")
# # write.table(temp, "full_control_integrate_nofilt.txt", sep="\t", quote=F)
# # write.table(temp2, "ICGScell_control_integrate_nofilt.txt", sep="\t", quote=F)
# 
# #make T cell plot
# resolution=0.3
# gene_list=c("CD4", "CD8A")
# pdf(file = paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_", resolution, "_FLT3/UMAP_Seurat_Tcell.pdf"), width = 10, height = 5)
# par(mar=c(1, 1, 1, 1))
# layout(matrix(1:2, ncol=2, nrow=1, byrow=TRUE))
# 
# for(n in 1:length(gene_list)){
#   gene=gene_list[n]
#   #filt_exp=immune.combined2@assays[["RNA"]]@data[,Seu_ICGS[,4]]
#   filt_exp=immune.combined2@assays[["RNA"]]@data
#   ind=which(row.names(filt_exp)==gene)
#   if(length(ind)!=0){
#     df=data.frame(cell=Seu_ICGS[,1], exp=filt_exp[ind,])
#     rbPal <- colorRampPalette(c('grey','red'))
#     df$Col <- rbPal(10)[as.numeric(cut(df$exp,breaks = 10))]
#     plot(Seu_ICGS_plot, pch = 16, cex=0.2, col = df$Col, main = gene, axes=FALSE, frame.plot=TRUE)
#     Axis(side=1, labels=FALSE, tick=F)
#     Axis(side=2, labels=FALSE, tick=F)
#   }else{
#     print(paste0("Gene ", gene, " not found"))
#   }
# }
# dev.off()
# 
# #make DC plot
# resolution=0.3
# gene_list=c("IL3RA", "CLEC4C", "NRP1", "FCER1A")
# pdf(file = paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_", resolution, "_FLT3/UMAP_Seurat_DC.pdf"), width = 6, height = 6)
# par(mar=c(1, 1, 1, 1))
# layout(matrix(1:4, ncol=2, nrow=2, byrow=TRUE))
# 
# for(n in 1:length(gene_list)){
#   gene=gene_list[n]
#   #filt_exp=immune.combined2@assays[["RNA"]]@data[,Seu_ICGS[,4]]
#   filt_exp=immune.combined2@assays[["RNA"]]@data
#   ind=which(row.names(filt_exp)==gene)
#   if(length(ind)!=0){
#     df=data.frame(cell=Seu_ICGS[,1], exp=filt_exp[ind,])
#     rbPal <- colorRampPalette(c('grey','red'))
#     df$Col <- rbPal(10)[as.numeric(cut(df$exp,breaks = 10))]
#     plot(Seu_ICGS_plot, pch = 16, cex=0.2, col = df$Col, main = gene, axes=FALSE, frame.plot=TRUE)
#     Axis(side=1, labels=FALSE, tick=F)
#     Axis(side=2, labels=FALSE, tick=F)
#   }else{
#     print(paste0("Gene ", gene, " not found"))
#   }
# }
# dev.off()
# 
# #Griffin pDC
# Griffin=c("JCHAIN",
# "IRF8",
# "TCF4",
# "APP",
# "CCDC50",
# "UGCG",
# "IRF4",
# "ITM2C",
# "GZMB",
# "MAPKAPK2",
# "LILRA4",
# "IRF7",
# "SERPINF1",
# "PLD4",
# "TGFBI",
# "SPIB",
# "SLC15A4",
# "LDLRAD4",
# "SMPD3",
# "MPEG1",
# "PPP1R14B",
# "PLAC8",
# "NIBAN3",
# "SCT",
# "CSF2RB",
# "TPM2",
# "BCL11A",
# "LGMN",
# "MAP1A",
# "SEL1L3",
# "CD74",
# "RRBP1",
# "HSP90B1",
# "SOX4",
# "RUNX2",
# "P2RY14",
# "CXXC5",
# "RGS1",
# "CTSB",
# "C12orf75",
# "FLNB",
# "UBE2J1",
# "RNASE6",
# "ZDHHC17",
# "CLEC4C",
# "PIK3AP1",
# "CD2AP",
# "RASD1",
# "MZB1",
# "TXNDC5")
# 
# #HummanImmuneCells
# HumanImmune=c("RNASET2",
# "DAPK2",
# "CYP46A1",
# "GPM6B",
# "LTK",
# "ZFAT",
# "DRD4",
# "SCT",
# "SLC12A3",
# "SIDT1",
# "TRAF4",
# "TTC39A",
# "DNAJA1",
# "RGS1",
# "LEPREL1",
# "GZMB",
# "SNTA1",
# "ASIP",
# "SMPD3",
# "CRYM",
# "PTPRS",
# "CADM4",
# "TSPAN13",
# "AEBP1",
# "PTGDS",
# "PALD1",
# "MAP2K6",
# "UNC93B1",
# "CUX2",
# "PTK7",
# "SLC4A3",
# "PROC",
# "NPC2",
# "BCL11A",
# "TGFBI",
# "PDZRN3",
# "CAT",
# "FAM213A",
# "C12orf44",
# "PFKFB2",
# "TOX2",
# "KCNK17",
# "LAMP5",
# "IGFLR1",
# "KCNA5",
# "PLVAP",
# "NAPSA",
# "SERPINF1",
# "ALOX5AP",
# "AMIGO2",
# "SRP14",
# "TSPAN3",
# "IRF8",
# "NPC1",
# "EIF4A3",
# "APP",
# "EPHA2",
# "ADC",
# "SDE2",
# "PLAC8",
# "RPS3A",
# "TNFRSF21",
# "UGCG",
# "C9orf142",
# "KIRREL3",
# "NGLY1",
# "CCDC50",
# "EPHB1",
# "PPM1J",
# "EEF1A1",
# "MX1",
# "WDR66",
# "NEK8",
# "CYB561A3",
# "RPS6KA4",
# "IDH3A",
# "PLD4",
# "MAP1A",
# "PHB",
# "TMIGD2",
# "SERPINF2",
# "KCTD5",
# "LDLRAD4",
# "ABHD15",
# "GPR183",
# "CLIC3",
# "MCC",
# "LRRC34",
# "PPP1R14B",
# "VEGFB",
# "P2RY14",
# "ZDHHC14",
# "GAPT",
# "FUT7",
# "HIGD1A",
# "CLN8",
# "SPNS3",
# "GAS6",
# "SLC35F3",
# "LRRC26",
# "CA13",
# "IL3RA",
# "IFNLR1",
# "IRF7",
# "KRT5",
# "S100A3",
# "C16orf93",
# "TCF4",
# "MPEG1",
# "CLEC4C",
# "TPM2",
# "SCAMP5",
# "ADAT3",
# "LCNL1",
# "PLXNA4",
# "AC007381.2",
# "AC011893.3",
# "SCAMP4",
# "RP11-71G12.1",
# "RP11-117D22.2",
# "LINC00865",
# "TLR9",
# "HHIP-AS1",
# "RP11-73G16.2",
# "AP000783.1",
# "RP11-783K16.5",
# "RP4-647C14.2",
# "SMIM6",
# "RP1-302G2.5",
# "BLNK",
# "C1orf186",
# "MAPKAPK2",
# "PLS3",
# "CD4",
# "RP11-542M13.3",
# "PLP2",
# "SPIB",
# "LILRA4",
# "NLRP7")
# 
# 
# # apply new cells to UMAP
# setwd("~/Documents/Projects/Data_files_temp/Donor")
# immune.query=readRDS("PB01_viable/Seurat.rds")
# setwd("~/Documents/Projects/Data_files_temp/Donor")
# immune.combined=readRDS("Seurat_Integration_0.5/Seurat_int_0.5.rds")
# 
# # immune.query <- MapQuery(anchorset = immune.anchors, reference = immune.combined, query = Q1, 
# #                            refdata = list(orig.ident = "orig.ident"), reference.reduction = "pca", reduction.model = "umap")
# 
# immune.anchors2 <- FindTransferAnchors(reference = immune.combined, query = immune.query, 
#                                         dims = 1:30, reference.reduction = "pca")
# 
# immune.query <- TransferData(anchorset = immune.anchors, reference = immune.combined, query = immune.query, 
#                                refdata = list(orig.ident = "orig.ident"))
# immune.query <- IntegrateEmbeddings(anchorset = immune.anchors, reference = immune.combined, 
#                                       query = immune.query, new.reduction.name = "ref.pca")
# immune.query <- ProjectUMAP(query = immune.query, query.reduction = "ref.pca", reference = immune.combined, 
#                               reference.reduction = "pca", reduction.model = "umap")
