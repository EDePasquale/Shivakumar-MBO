rm(list=ls())

#load and organized the file
library(readxl)
dat_pre <- read_xlsx("STM_ somalogic table_2017.xlsx")
table(dat_pre[,1])
colnames(dat_pre)[1]="Diagnosis"
dat <- dat_pre[dat_pre$Diagnosis%in%c("BA","NC"),]
dat$Group <- c(
  rep("BA",175),
  rep("NC",9)
)
dat$Group <- factor(dat$Group, levels=unique(dat$Group))

#calculate mean protein level by group
library(dplyr)

Meanbygroup <- apply(dat[,-c(1,1129)], 2, function(x) tapply(x,dat$Group, mean))
Meanbygroup <- t(Meanbygroup)
Meanbygroup[1:4,]
Meanbygroup <- as.data.frame(Meanbygroup)
str(Meanbygroup)

#calculate log2FC
log2FC <- Meanbygroup %>% 
  mutate(`BAvsNC` = log2(BA/NC))

#perform Student’s t test
testFunction <- function (x) {
  return(tryCatch(t.test(x ~ dat$Group,var.equal=TRUE), error=function(e) NULL))
}

p <- lapply(dat[,-c(1,1129)], testFunction)
v1 <- sapply(p, function(x) x$p.value)
names(v1) <-  names(p)
y <-lapply(v1, function(i)replace(i, length(i) == 0, 1))
tt <-do.call(rbind, unname(Map(cbind, Names = names(y), y)))
tt<- as.data.frame(tt)
tt[,2] <- as.numeric(tt[,2])
colnames(tt)[2] <- "pvalue"

#calculate FDR
tt$padj = p.adjust(tt$pvalue, method ="BH", n=1127)

#save the result
Result = cbind(log2FC, tt[,-1])
colnames(Result)[4]="log2FC"
da <- read_xlsx("target_GeneSymbol_list.xlsx")
gene_list <- t(da)
gene_list <- cbind(rownames(gene_list),gene_list)
colnames(gene_list) <- gene_list[2,]
gene_list <- as.data.frame(gene_list[-(1:2),])
Full_result <- cbind(gene_list,Result[,-1])
library(openxlsx)
library(xlsx)
wb = createWorkbook()
sheet = createSheet(wb)
addDataFrame(Full_result, sheet=sheet, startColumn=1, row.names=FALSE)
saveWorkbook(wb, "Full_result.xlsx")
