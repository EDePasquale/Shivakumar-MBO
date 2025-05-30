# Load libraries
library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

# Load Seurat object and set working directory
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5/")
M<-readRDS("Seurat_8_MBO_withMeta_names.rds")
M <- SetIdent(M, value = M@meta.data$cluster_names)

# Add cluster names (change transitional to intermediate)
cluster_names<-c("Mesenchymal",  
                 "Epithelial",   
                 "Epithelial",   
                 "Endothelial",  
                 "Epithelial",   
                 "Intermediate 1",
                 "Intermediate 2",
                 "Intermediate 3",
                 "Mesenchymal",  
                 "Mesenchymal",  
                 "Intermediate 4", 
                 "Intermediate 5", 
                 "Epithelial",   
                 "Other",          
                 "Epithelial"  )
M@meta.data[["merged_clusters"]] <- mapvalues(M@meta.data[["cluster_names"]], from=levels(M@meta.data$cluster_names), to=cluster_names)
Idents(M) = M$merged_clusters

pdf(file = "UMAP_All_Names_Merged.pdf", width = 10, height = 9)
par(mar=c(2, 2, 2, 2))
DimPlot(M, reduction="umap", group.by = "merged_clusters", label=T, raster=F, repel=T)
dev.off()

pdf(file = "UMAP_Split_Names_Merged.pdf", width = 17, height = 9)
par(mar=c(2, 2, 2, 2))
  DimPlot(M, reduction="umap", group.by = "merged_clusters", label=F, raster=F, repel=T, split.by = "Diagnosis")
dev.off()

# Make proportion plots by diagnosis
data_long=as.data.frame(cbind(Diagnosis=M@meta.data[["Diagnosis"]], Cluster=as.character(M@meta.data[["merged_clusters"]])))
data=as.data.frame.matrix(table(data_long))

# Read and wrangle data
data=t(data)
data <- sweep(data, 2, colSums(data), "/")*100
write.table(data, "Cell_Type_Proportions.txt", sep="\t", quote=F)

gg_color_hue <- function(n) { # ggplot colors
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
my_colors=gg_color_hue(9)

myColors=as.data.frame(cbind(Cluster_Name=levels(M@meta.data[["merged_clusters"]]), color=my_colors))
myColors=myColors[order(myColors$Cluster_Name, decreasing=F),]
my_colors=myColors$color

# Make plot split like cluster (see Hiro example)
my_colors_fade=alpha(my_colors, 0.2)
idx <- order(c(seq_along(my_colors_fade), seq_along(my_colors)))
new_colors=unlist(c(my_colors_fade, my_colors))[idx]

pdf("CellTypeFrequencies_byDiagnosis_Hiro_merged.pdf", width = 6, height = 9)
par(mar = c(6,10,2,2), xpd = T)
  barplot(height=t(as.matrix(data[nrow(data):1,])), beside=T, col = rev(new_colors), horiz=T, las=1, xlab="Proportion (%)")
  legend(x="bottomright", legend=c("Control", "Biliary Atresia"), fill=new_colors[1:2])
dev.off()

# Make it with error bars
data_sample_long=as.data.frame(cbind(Donor=M@meta.data[["Donor"]], Cluster=as.character(M@meta.data[["merged_clusters"]])))
data_sample=as.data.frame.matrix(table(data_sample_long))
data_sample=data_sample[c(1,4,5,6,2,3,7,8),]
data_sample=t(data_sample)
data_sample <- sweep(data_sample, 2, colSums(data_sample), "/")*100
data_sample=t(data_sample)
NC_SD=apply(data_sample[1:4,], 2, sd)
BA_SD=apply(data_sample[5:8,], 2, sd)
data_SD=t(as.matrix(cbind(BA_SD, NC_SD)))

#A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(y+upper, x,y, x,angle=90, code=3, length=length, ...)
}

pdf("CellTypeFrequencies_byDiagnosis_Hiro_merged_errorbars.pdf", width = 6, height = 9)
par(mar = c(6,10,2,4), xpd = T)
  X=barplot(height=t(as.matrix(data[nrow(data):1,])), beside=T, col = rev(new_colors), horiz=T, las=1, xlab="Proportion (%)")
  error.bar(X,t(as.matrix(data[nrow(data):1,])), data_SD)
  legend(x="bottomright", legend=c("NC", "BA"), fill=c("gray78", "gray32"))
dev.off()

##### 
# Do it with ggplot
data_sample=rbind(data_sample[5:8,], data_sample[1:4,])
data_sample=cbind(sample=row.names(data_sample), data_sample)
new_data_sample=reshape2::melt(data_sample, id.vars = "sample")
colnames(new_data_sample)=c("Sample", "Cluster", "Frequency")
new_data_sample=new_data_sample[-c(1:8),]
new_data_sample2=cbind(new_data_sample, Group=rep(c(rep("BA", 4), rep("NC", 4)),9))
new_data_sample3=new_data_sample2[,c(2:4)]
row.names(new_data_sample3)<-NULL
new_data_sample3$Frequency<-as.numeric(new_data_sample3$Frequency)
knitr::kable(head(new_data_sample3))
#   |Cluster     | Frequency|Group |
#   |:-----------|---------:|:-----|
#   |Endothelial | 16.080247|BA    |
#   |Endothelial | 17.309472|BA    |
#   |Endothelial |  5.285776|BA    |
#   |Endothelial | 14.035793|BA    |
#   |Endothelial |  2.107728|NC    |
#   |Endothelial |  6.536388|NC    |

new_colors_2=new_colors[c(2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17)]

p1=ggplot(new_data_sample3, aes(x=reorder(Cluster, Frequency), y=Frequency, fill=interaction(Group, Cluster))) + 
  geom_boxplot(outlier.color=NA) +
  geom_dotplot(binaxis="y", binwidth=0.5, stackdir="center", dotsize=2, position = position_dodge(width = 0.75)) +
  coord_flip() +
  scale_fill_manual(values=new_colors_2) +
  xlab(NULL) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  +
  theme(legend.position="none")

pdf("CellTypeFrequencies_byDiagnosis_Hiro_merged_ggplot.pdf", width = 6, height = 5)
par(mar = c(6,10,2,4), xpd = T)
print(p1)
dev.off()

#######

# Make proportion plots by diagnosis
data_long=as.data.frame(cbind(Donor=M@meta.data[["Donor"]], Cluster=as.character(M@meta.data[["merged_clusters"]])))
data=as.data.frame.matrix(table(data_long))

data=data[c(1,4,5,6,2,3,7,8),]

#Intermediate 4 (column 6)
NC=data[c(1:4),6]
BA=data[c(5:8),6]

data=t(data)
data <- sweep(data, 2, colSums(data), "/")*100
data=t(data)

for(i in 1:ncol(data)){
  NC=data[c(1:4),i]
  BA=data[c(5:8),i]
  x=t.test(NC,BA)
  print(paste0(colnames(data)[i], ": ", x$p.value))
}

# [1] "Endothelial: 0.156305151301009"
# [1] "Epithelial: 0.0404582534052831"
# [1] "Intermediate 1: 0.232326013679259"
# [1] "Intermediate 2: 0.00216455708459744"
# [1] "Intermediate 3: 0.0547328569571133"
# [1] "Intermediate 4: 0.0268482859002916"
# [1] "Intermediate 5: 0.223773342736508"
# [1] "Mesenchymal: 0.432762613011513"
# [1] "Other: 0.675707336841021"
