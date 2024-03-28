# library(devtools)
# install_github("wjawaid/enrichR")
library(enrichR)
library(Seurat)
library(ggplot2)
library(scales)
library(plyr)
library(dplyr)

websiteLive <- getOption("enrichR.live") #https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

if (websiteLive) dbs <- listEnrichrDbs() # all available gene sets

setwd("/Volumes/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5/")
M<-readRDS("Seurat_8_MBO_withMeta_names.rds")
M <- SetIdent(M, value = M@meta.data$cluster_names)

# Add cluster names
cluster_names<-c("Mesenchymal",  
                 "Epithelial",   
                 "Epithelial",   
                 "Endothelial",  
                 "Epithelial",   
                 "Intermediate 1", # MEG3 Transitional
                 "Intermediate 2", # LGALS1 Transitional
                 "Intermediate 3", # CYP3A5 Transitional
                 "Mesenchymal",  
                 "Mesenchymal",  
                 "Intermediate 4", # DCN Transitional
                 "Intermediate 5", # TFF1 Transitional
                 "Epithelial",   
                 "Other",          
                 "Epithelial")
M@meta.data[["merged_clusters"]] <- mapvalues(M@meta.data[["cluster_names"]], from=levels(M@meta.data$cluster_names), to=cluster_names)
Idents(M) = M$merged_clusters

M<-subset(M, idents=c("Intermediate 1", "Intermediate 2", "Intermediate 3", "Intermediate 4", "Intermediate 5"))
M.markers=FindAllMarkers(M, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 1, return.thresh = 0.05)

M.markers=M.markers[which(M.markers$p_val_adj <=0.05),]

dbs <- c("MSigDB_Hallmark_2020")
if (websiteLive) {
  for(i in 1:5){ # TODO: CANT RUN IN LOOP, MUST RUN MANUALLY!
    my_genes=M.markers$gene[which(M.markers$cluster==paste0("Intermediate ", i))]
    print(my_genes) #printing as different for each i (GOOD)
    x=enrichr(my_genes, dbs) #when done outside of a loop this works... (BAD)
    print(x) #showing same results with each i (BAD)
    assign(paste0("enriched_", i), x)
  }
}

# 1, 2, and 4 have EMT as top (as predicted)

## set the levels in order we want
my_order=enriched_1[["MSigDB_Hallmark_2020"]][,1]
enriched_1[["MSigDB_Hallmark_2020"]][,1]<-factor(enriched_1[[1]][,1], levels=my_order)
my_order=enriched_2[["MSigDB_Hallmark_2020"]][,1]
enriched_2[["MSigDB_Hallmark_2020"]][,1]<-factor(enriched_2[[1]][,1], levels=my_order)
my_order=enriched_3[["MSigDB_Hallmark_2020"]][,1]
enriched_3[["MSigDB_Hallmark_2020"]][,1]<-factor(enriched_3[[1]][,1], levels=my_order)
my_order=enriched_4[["MSigDB_Hallmark_2020"]][,1]
enriched_4[["MSigDB_Hallmark_2020"]][,1]<-factor(enriched_4[[1]][,1], levels=my_order)
my_order=enriched_5[["MSigDB_Hallmark_2020"]][,1]
enriched_5[["MSigDB_Hallmark_2020"]][,1]<-factor(enriched_5[[1]][,1], levels=my_order)

mycolor=hue_pal()(15) # original number of samples
graph_colors=mycolor[c(6,7,8,11,12)]

#ggplot(enriched_1[["MSigDB_Hallmark_2020"]][1:5,], aes(x=Term, y=P.value)) + geom_bar(fill=mycolor[6], stat="identity")


test=rbind(cbind(enriched_1[["MSigDB_Hallmark_2020"]][1:5,], Graph=rep("Intermediate 1", 5)),
           cbind(enriched_2[["MSigDB_Hallmark_2020"]][1:7,], Graph=rep("Intermediate 2", 7)),
           cbind(enriched_3[["MSigDB_Hallmark_2020"]][1:5,], Graph=rep("Intermediate 3", 5)),
           cbind(enriched_4[["MSigDB_Hallmark_2020"]][1:5,], Graph=rep("Intermediate 4", 5)),
           cbind(enriched_5[["MSigDB_Hallmark_2020"]][1:5,], Graph=rep("Intermediate 5", 5)))

# Create -log pvalue column!!
test=cbind(test, neglogP.value=(-log10(test$P.value)))


pdf(file = "Enrichr_neglogPvalue_0.05.pdf", width = 6, height = 3)
par(mar=c(10, 6, 4,4))
test %>%
  ggplot(aes(x=Term, y=neglogP.value, fill=Graph)) + 
  scale_fill_manual(values=graph_colors) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  facet_grid(~Graph, scales="free", space="free") + 
  theme_bw() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="none", axis.title.x = element_blank()) + 
  labs(y = "-log10 P-value")
dev.off()
