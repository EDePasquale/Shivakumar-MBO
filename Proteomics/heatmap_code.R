rm(list = ls())
#load selected gene list
load("List.Rdata")

#filter protein level data with selected genes
load("Protein level data.Rdata")

dat <- t(dat)
da <- dat[match(List,rownames(dat)),]
da2 <- matrix( 
  as.numeric(da), ncol = 184,
  dimnames = list(rownames(da),
                  colnames(da)))
str(da2)

group <- data.frame(group=c(
  rep("BA",175),
  rep("NC",9)
))
rownames(group)=colnames(da2)
group$group <- factor(group$group,levels = c("NC","BA"))
str(group)


# make and save the heatmap
library(pheatmap)

p <- pheatmap(log2(da2),scale = "row",
              show_colnames =F,show_rownames = T, 
              cluster_cols = T,
              cluster_rows = T,
              annotation_col=group,
              breaks = seq(-1.5,1.5,length.out = 100),
              main = "Heatmap")
pdf(file = "Heatmap_proteins_EMT.pdf",width = 8,height = 6)
print(p)
dev.off()


