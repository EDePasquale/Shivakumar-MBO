#install.packages('SoupX')
library(SoupX)
#BiocManager::install("DropletUtils")
library(DropletUtils)

path_root_1="/data/livsyscom/5_MBO_scRNAseq_Oct2021/cellranger6/"
path_root_2="/data/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/CellRanger/"

exp_list_1=c("10X-Bezerra-CO15-MBO-20210816-3v3-1hg",
             "10X-Bezerra-CO36-MBO-20211103-3v3-1hg",
             "10X-Bezerra-CO-53-MBO-NC-20211201-3v3-1hg",
             "10X-Bezerra-CO55-MBO-20211027-3v3-1hg",
             "10X-Bezerra-CO21_MBO-BA-20220621",
             "10X-Bezerra-CO-022-MB0-20211117-3v3-1hg",  
             "10X-Bezerra-CO-68-MBO-20220302-3v3",
             "10X-Bezerra-CO72-MBO-BA-p10-20220509")
exp_list_2=c("CO15-MBO-NC-20210816",
           "CO36-MBO-NC-20211103",
           "CO53-MBO-NC-20211201",
           "CO55-MBO-NC-20211027",
           "CO21_MBO-BA-20220621",
           "CO22-MBO-BA-20211117",
           "CO68-MBO-BA-20220302",
           "CO72-MBO-BA-20220509")


for(i in 1:length(exp_list_1)){

  final_path=paste0(path_root_1, exp_list_1[i])
  setwd(final_path) #folder that contains outs
  
  #Decontaminate one channel of 10X data mapped with cellranger by running:
  sc = load10X(paste0(final_path, "/outs"))
  
  setwd(paste0(path_root_2, exp_list_2[i]))
        
  pdf(file = paste0(final_path, ".pdf"), width = 11, height = 8.5)
  sc = autoEstCont(sc)
    print(sc)
  dev.off()
  
  # Save estimated rho
  temp=sc[["fit"]][["rhoEst"]]
  names(temp)=exp_list_2[i]
  write.table(temp, "est.rho.txt", sep="\t", quote=F, row.names = T, col.names = F)
  
  out = adjustCounts(sc)
  #out will then contain a corrected matrix to be used in place of the original table of counts in downstream analyses.

  cntSoggy = Matrix::rowSums(sc$toc > 0)
  cntStrained = Matrix::rowSums(out > 0)
  mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 100)
  write.table(mostZeroed, paste0(exp_list_2[i], "_mostZeroed.txt"), sep="\t")
  #mostly cell type specific genes
  
  leastZeroed = tail(sort(Matrix::rowSums(sc$toc > out)/Matrix::rowSums(sc$toc > 0)), n = 100)
  write.table(leastZeroed, paste0(exp_list_2[i], "_leastZeroed.txt"), sep="\t")
  #general genes
  
  DropletUtils:::write10xCounts("./strainedCounts", out, version="3")
  
}