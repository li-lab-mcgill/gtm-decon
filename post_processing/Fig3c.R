#Figure 3c: Deconvolution of GEO-Fadista dataset - using PI_Segerstolpe as training dataset

library(pheatmap)
library(reshape2)
library(ggplot2)
library(ggpubr)

library(forcats)
library(dplyr)
library(Rmisc)

f <- read.csv(file="../data/GEO_Fadista_5topics_metagene_normalized.csv",header=FALSE)
mat <- as.matrix(f[,c(1:9,12)])
colnames(mat)<-c("acinar","alpha","beta","delta","ductal","endothelial","epsilon","gamma","mast","PSC")
annot_row = data.frame(Disease=factor(c("N","na","N","na","na","N","na","N","na","na","T2D","na","N","N","N","N","IGT","N","N","na","N", "na", "N",
                                        "N","IGT","N", "IGT", "N",  "T2D", "N",  "N",  "T2D", "N", "N",  "IGT", "N", "N", "IGT",  "N",  "N",   "na",  "T2D",  "T2D",  "N",  "N",  "T2D",
                                        "N", "IGT",  "na", "N", "IGT",  "N", "IGT",  "IGT",  "T2D", "T2D",  "N",  "N",  "N",  "N",  "IGT",  "N",  "N",  "N",  "IGT",  "N",  "N",  "N",  "N",
                                        "N", "N",  "N",  "N",  "N",  "N",  "IGT",  "IGT",  "T2D",  "na",  "N",  "N",  "T2D",  "IGT",  "T2D",  "N",  "N",  "IGT",  "N",  "N")))

#rownames(mat) <- paste0("row_", seq(nrow(f)))
rownames(annot_row) <- rownames(f)
ann_colors = list(Disease=c(N = "yellow", IGT = "lightgreen", T2D = "darkgreen", na = "white" ),
                  HbA1c = c("orange","blue"))


matrix <-t(apply(mat,1, function(x) x/sum(x))) #Change each element by value / rowSum
rownames(matrix) <- rownames(f)

png("../data/Fig3c.png", units="in", width=7, height=7, res=300)
pheatmap(matrix, annotation_row = annot_row, annotation_colors = ann_colors, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,cluster_cols=F,cluster_rows=T,fontsize = 14, border_color = NA) 
dev.off()