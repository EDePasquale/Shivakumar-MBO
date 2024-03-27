# Libraries
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(chorddiag)  #devtools::install_github("mattflor/chorddiag")
library(scales)

TGFB=c("CPI-CC0104F2A96", "CPI-CC072007F09", "CPI-CS0533A30DC", "CPI-CS0C159B54D", "CPI-CS02899BF8B", "CPI-CS085EF890A", "CPI-SC098C06330", "CPI-SC0984894BC")
WNT=c("CPI-SS0B2D5BA9", "CPI-SS0D81C6C1D", "CPI-SS01C2BA5C8", "CPI-SS0E32BEF02", "CPI-SS002673EFE", "CPI-SS0B8083B44", "CPI-SS082AD7B6A")
NOTCH=c("CPI-SS0C4D54E39", "CPI-SS0B87FE6E0", "CPI-SS060B82BBC", "CPI-SS0B99F22A2")

####################
#                  #
# All Interactions #
#                  #
####################

# Load dataset
my_data <- read.table("/Volumes/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5/CellPhoneDB/All/Shivakumar_DEG_results/relevant_interactions_edit_names.txt", sep="\t", header=F)
my_data=my_data[c(1,which(my_data[,1] %in% TGFB)),] # REDUCTION to TGFB interaction set
my_data2 <- my_data[,-c(1:11)]
names_list=apply(my_data2[1,], 1, strsplit, "_")
my_data3=as.data.frame(do.call(rbind, names_list[[1]]))
colnames(my_data3)<-c("rowname", "key")
a=as.matrix(my_data2[2:nrow(my_data2),])
a<- matrix(as.numeric(a),    # Convert to numeric matrix
                    ncol = ncol(a))
my_data3=cbind(my_data3, value=colSums(a))

# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# color palette
mycolor=hue_pal()(9)
mycolor=mycolor[c(6,7,3,2,5,4,1,9,8)]

pdf(file = "/Volumes/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5/CellPhoneDB/All/Shivakumar_DEG_results/Chord_allInteractions_TGFB.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
# Base plot
chordDiagram(
  x=my_data3,
  grid.col = mycolor,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE)

# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 3.2, 
      labels = sector.index, 
      facing = "bending", 
      cex = 0.8
    )
    
    # Add graduation on axis
    circos.axis(
      h = "top", 
      #major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)), 
      #minor.ticks = 5, 
      #major.tick.percentage = 0.5,
      labels.niceFacing = FALSE)
  }
)
dev.off()

#############################
#                           #
# Non-Integrin Interactions #
#                           #
#############################

# Load dataset
my_data <- read.table("/Volumes/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5/CellPhoneDB/All/Shivakumar_DEG_results/relevant_interactions_edit_names.txt", sep="\t", header=F)
my_data=my_data[c(1,which(my_data[,1] %in% TGFB)),] # REDUCTION
my_data=my_data[c(1,which(my_data$V11=="FALSE")),]
my_data2 <- my_data[,-c(1:11)]
names_list=apply(my_data2[1,], 1, strsplit, "_")
my_data3=as.data.frame(do.call(rbind, names_list[[1]])) #not sure why this would be rbind and not cbind...
colnames(my_data3)<-c("rowname", "key")
a=as.matrix(my_data2[2:nrow(my_data2),])
a<- matrix(as.numeric(a),    # Convert to numeric matrix
           ncol = ncol(a))
my_data3=cbind(my_data3, value=colSums(a))

# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# color palette
mycolor=hue_pal()(9)
mycolor=mycolor[c(6,7,3,2,5,4,1,9,8)]

pdf(file = "/Volumes/GI-Informatics/DePasquale/Projects/Shivakumar_MBO/Seurat_Integration_0.5/CellPhoneDB/All/Shivakumar_DEG_results/Chord_nonIntegrinInteractions_TGFB.pdf", width = 8, height = 8)
par(mar=c(2, 2, 2, 2))
# Base plot
chordDiagram(
  x=my_data3,
  grid.col = mycolor,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE)

# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 3.2, 
      labels = sector.index, 
      facing = "bending", 
      cex = 0.8
    )
    
    # Add graduation on axis
    circos.axis(
      h = "top", 
      #major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)), 
      #minor.ticks = 5, 
      #major.tick.percentage = 0.5,
      labels.niceFacing = FALSE)
  }
)
dev.off()

