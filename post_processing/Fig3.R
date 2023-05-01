## Figure 3a: Visualing top20 genes of genes-by-CTS-topics matrix, generated using 5topics per cell type, on PI_Segerstolpe scRNA-seq dataset
library(pheatmap)
library(reshape2)
library(ggplot2)
library(ggpubr)

library(forcats)
library(dplyr)
library(Rmisc)

# Get genesets for pancreatic tissue from CellMarker Database
cellmarker <- read.csv(file="../data/Human_cell_markers_geneset.txt",sep="#")
cellMarker_gs <- split(x = cellmarker$GeneSymbol, f = cellmarker$TissueCellName)

f <-read.csv(file="../data/trainData_5topics_phi_normalized.csv",header=FALSE)
genes <-read.csv(file="../data/genes.txt",header=FALSE)

df <- cbind(genes,f)
#Drop the first second and third column - they dont correspond to phi values
df<-df[-(2:3)]
#Get top20 genes per column
datalist = list()
topn_df <- data.frame()
for(i in 2:ncol(df)) {       # for-loop over columns
  subset_df <- df %>% arrange(desc(df[,i])) %>% slice(1:20)
  datalist[[i]]<-subset_df
}
topn_df <- do.call(rbind,datalist)
max_val <- max(select(topn_df,2:ncol(topn_df)))
topn_df_rescaled <- topn_df
topn_df_rescaled[2:ncol(topn_df)] <- topn_df[2:ncol(topn_df)]/max_val

# For 5 topics
matrix <- as.data.frame(select(topn_df,2:46,57:61)) #Column selection - removing unncessary cell types
matrix <- matrix[-c(901:1100,1201:1400),]  #Row selection - removing unnecessary cell types
max_val <- max(matrix)
matrix <- matrix/max_val

# Getting annotation list for all top20 genes in the sets
a <-topn_df[1:100,]$V1
annot_list <- a %in% cellMarker_gs[["Pancreatic islet - Acinar cell"]]
a <-topn_df[101:200,]$V1
b <- a %in% cellMarker_gs[["Pancreatic islet - Alpha cell"]]
annot_list <- append(annot_list,b)
a <-topn_df[201:300,]$V1
b <- a %in% cellMarker_gs[["Pancreatic islet - Beta cell"]]
annot_list <- append(annot_list,b)
a <-topn_df[301:400,]$V1
b <- a %in% cellMarker_gs[["Pancreas - Delta cell"]]
annot_list <- append(annot_list,b)
a <-topn_df[401:500,]$V1 #Only 1 annotation available for ductal cell
#b <- c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA")
b <- a %in% cellMarker_gs[["Pancreas - Ductal cell"]]
annot_list <- append(annot_list,b)
a <-topn_df[501:600,]$V1 # no annotations available for endothelial cell
b <- c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",
       "NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",
       "NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",
       "NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",
       "NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA")
annot_list <- append(annot_list,b)
a <-topn_df[601:700,]$V1 # no annotations available for epsilon cell
annot_list <- append(annot_list,b)
a <-topn_df[701:800,]$V1 # no annotations available gamma cell
annot_list <- append(annot_list,b)
a <-topn_df[801:900,]$V1 # no annotations available mast cell
annot_list <- append(annot_list,b)
a <-topn_df[1101:1200,]$V1 # no annotations available PSC cell
annot_list <- append(annot_list,b)

#Acinar from PanglaoDB
m_list <- c("PRSS1","KLK1","PNLIP","CTRC","ALDOB","REG3A","SERPINA3","PRSS3","REG1B","CFB","GDF15","MUC1","C15ORF48","DUOXA2","AKR1C3","OLFM4","GSTA1","LGALS2","PDZK1IP1","RARRES2","CXCL17","UBD","GSTA2","ANPEP","LYZ","ANGPTL4","CTRB1","RBPJL","PTF1A","CELA3A","SPINK1","ZG16","CEL","CELA2A","CPB1","CELA1","PNLIPRP1","RNASE1","AMY2B","CPA2","CPA1","CELA3B","CTRB2","PLA2G1B","PRSS2","CLPS","REG1A","SYCN")
a <-topn_df[1:100,]$V1
annot_list_panglaodb <- a %in% m_list

a <-topn_df[101:200,]$V1
#Alpha from PanglaoDB
m_list <- c("FXYD5","LDB2","MAFB","PCSK2","CHGA","FEV","SCGB2A1","FAP","DPP4","GPR119","PAX6","NEUROD1","LOXL4","PLCE1","GC","KLHL41","SLC30A8","PTGER3","RFX6","SMARCA1","PGR","IRX1","UCP2","RGS4","GLS","KCNK16","TTR","GLP1R","ARX","POU3F4","NKX2-2","RESP18","PYY","SLC38A5","TM4SF4","CRYBA2","SH3GL2","PCSK1","PRRG2","IRX2","ALDH1A1","PEMT","SMIM24","F10","SCGN","GCG","NKX6-1")
b <- a %in% m_list
annot_list_panglaodb <- append(annot_list_panglaodb,b)

a <-topn_df[201:300,]$V1
#Beta from PanglaoDB
m_list <- c("NKX6-1","NKX2-2","MAFA","INS","SLC2A2","PDX1","GJD2","SH3GL2","HOPX","ADCYAP1","SCGB2A1","EDARADD","CASR","PFKFB2","MAFB","PAX6","RGS16","NPTX2","NEUROD1","ISL1","TGFBR3","SMAD9","SIX3","BMP5","PIR","STXBP5","DLK1","MEG3","GCGR","LMX1A","JPH3","CD40","HAMP","EZH1","NTRK1","NKX6-2","FXYD2","NPY","RIMS1","EFNA5","FFAR2","PAX4","IAPP","PCSK2","G6PC2","SLC30A8","PCSK1","SCGN","IGF2","SYT13","SIX2")
b <- a %in% m_list
annot_list_panglaodb <- append(annot_list_panglaodb,b)

a <-topn_df[301:400,]$V1
#Delta from PanglaoDB
m_list <- c("SST","FRZB","MS4A8","BAIAP3","CASR","BCHE","GABRB3","UNC5B","EDN3","PRG4","GHSR","PCSK1","GABRG2","POU3F1","BHLHE41","EHF","LCORL","ETV1","PDX1","LEPR","UCP2","GPC5-AS1","NPTX2","FXYD2","IAPP","KCNK16","SCGN","ISL1","HHEX","RESP18","PAX4","RBP4","PCSK9","FFAR4")
b <- a %in% m_list
annot_list_panglaodb <- append(annot_list_panglaodb,b)

a <-topn_df[401:500,]$V1 
#Ductal from PanglaoDB
m_list <- c("TFF1","CLDN1","PIGR","LGALS4","PERP","PDLIM3","WFDC2","SLC3A1","AQP1","ALDH1A3","VTCN1","KRT19","TFF2","KRT7","CLDN4","LAMB3","TACSTD2","CCL2","DCDC2","CXCL2","CTSH","S100A10","CFB","HNF1B","KRT20","MUC1","ONECUT1","AMBP","HHEX","ANXA4","SPP1","PDX1","SERPINA3","CFTR","GDF15","AKR1C3","MMP7","DEFB1","TSPAN8","CLDN10","SLPI","SERPINA5","SERPING1")
b <- a %in% m_list
annot_list_panglaodb <- append(annot_list_panglaodb,b)

a <-topn_df[501:600,]$V1 # no annotations available for endothelial cell
b <- c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",
       "NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA")
annot_list_panglaodb <- append(annot_list_panglaodb,b)

a <-topn_df[601:700,]$V1 
#Epsilon from PanglaoDB
m_list <- c("GHRL","ARX","HRH2","CALCR","SLC6A16","PCSK6","ADAMTS6","COL22A1","FAM124A","COL12A1","CD109","THSD4","CORIN","NKX2-2","ACSL1","FRZB","VTN","APOH","SPINK1","VSTM2L","SPTSSB","S100A6","KRT8","ANXA13","PHGR1","BMP4","HMGCS2","TM4SF5","OLFML3","ASGR1","COX6A2","NPY1R","FAXDC2","SLC7A9","MYO1A","C2ORF54","GHRLOS","BHMT","OPRK1","PTGER4","ITGA4","EYA4","XYLB","ELOVL2","AFF3")
b <- a %in% m_list
annot_list_panglaodb <- append(annot_list_panglaodb,b)

a <-topn_df[701:800,]$V1 
#Gamma (PP) from PanglaoDB
m_list <- c("PPY","ABCC9","FGB","ZNF503","MEIS1","LMO3","EGR3","CHN2","PTGFR","ENTPD2","AQP3","THSD7A","CARTPT","ISL1","PAX6","NEUROD1","APOBEC2","SEMA3E","SLITRK6","SERTM1","PXK","PPY2P","ETV1","ARX","CMTM8","GPC5-AS1","SCGB2A1","FXYD2","SCGN")
b <- a %in% m_list
annot_list_panglaodb <- append(annot_list_panglaodb,b)

a <-topn_df[801:900,]$V1 # no annotations available mast cell
b <- c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",
       "NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA")
annot_list_panglaodb <- append(annot_list_panglaodb,b)

a <-topn_df[1101:1200,]$V1 # no annotations available PSC cell
m_list <- c("COL6A1","THY1","TIMP1","COL6A3","SFRP2","COL1A2","TNFAIP6","TIMP3","COL1A1","SPARC","COL3A1","MGP","COL6A2","COL4A1","FN1","MMP11","TGFB1","INHBA","PDGFRA","NDUFA4L2","MMP14","CTGF","CYGB","KRT10","PDGFRB","DYNLT1","GEM","SPON2","RGS5")
b <- a %in% m_list
annot_list_panglaodb <- append(annot_list_panglaodb,b)


annot_row = data.frame(CellMarkerDB = factor(annot_list),
                       PanglaoDB = factor(annot_list_panglaodb))
ann_colors = list(CellMarkerDB=c("TRUE" = "red", "FALSE" = "white", "NA" = "lightgrey"),
                  PanglaoDB = c("TRUE" = "orange", "FALSE" = "white", "NA" = "lightgrey"))

rownames(annot_row)=rownames(matrix)
celltypes <- c("acinar","alpha","beta","delta","ductal","endothelial","epsilon","gamma","mast","PSC")
png("../data/Fig3a.png", units="in", width=7, height=7, res=300)
pheatmap(matrix,annotation_row = annot_row, annotation_colors = ann_colors, cluster_cols=F,cluster_rows=F,color=colorRampPalette(c("white","blue"))(50),show_rownames = F,labels_col = celltypes,fontsize = 12)
dev.off()


# # b) GSEA for the main cell types
# df <- cbind(genes,f)
# rownames(df)<-df$V1
# 
# #Acinar
# gseaDat <- df
# ranks <- gseaDat$V3
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
# 
# gseaDat <- df
# ranks <- gseaDat$V4
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
# 
# gseaDat <- df
# ranks <- gseaDat$V5
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
# 
# gseaDat <- df
# ranks <- gseaDat$V6
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
# 
# gseaDat <- df
# ranks <- gseaDat$V7
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
# 
# 
# #Alpha
# ranks <- gseaDat$V8
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# #head(fgseaRes[order(padj,-abs(NES)), ],2)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
# 
# ranks <- gseaDat$V9
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# #head(fgseaRes[order(padj,-abs(NES)), ],2)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
# 
# ranks <- gseaDat$V10
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# #head(fgseaRes[order(padj,-abs(NES)), ],2)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
# 
# ranks <- gseaDat$V11
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# #head(fgseaRes[order(padj,-abs(NES)), ],2)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
# 
# ranks <- gseaDat$V12
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# #head(fgseaRes[order(padj,-abs(NES)), ],2)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
# 
# 
# #Beta
# ranks <- gseaDat$V13
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# #head(fgseaRes[order(-abs(NES)), ],2)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
# 
# ranks <- gseaDat$V14
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# #head(fgseaRes[order(padj,-abs(NES)), ],2)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
# 
# ranks <- gseaDat$V15
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# #head(fgseaRes[order(padj,-abs(NES)), ],2)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
# 
# ranks <- gseaDat$V16
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# #head(fgseaRes[order(padj,-abs(NES)), ],2)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
# 
# ranks <- gseaDat$V17
# names(ranks) <- rownames(gseaDat)
# fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
# head(fgseaRes[order(pval), ],5)
# #head(fgseaRes[order(padj,-abs(NES)), ],2)
# plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
#                ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20))


# c,d,e,f) Deconvolution of GEO-Fadista dataset - using PI_Segerstolpe and PI_Baron as training datasets

f <- read.csv(file="../data/Fig3_GEO_Fadista_5topics_metagene_normalized.csv",header=FALSE)

mat <- as.matrix(f[,c(1:9,12)])
segerstolpe_mat <- mat
#mat <- as.matrix(f[,c(1:45,56:60)])
colnames(mat)<-c("acinar","alpha","beta","delta","ductal","endothelial","epsilon","gamma","mast","PSC")
#colnames(mat)<-c("acinar","2","3","4","5","alpha","2","3","4","5","beta","2","3","4","5","delta","2","3","4","5","ductal","2","3","4","5","endothelial","2","3","4","5","epsilon","2","3","4","5","gamma","2","3","4","5","mast","2","3","4","5","PSC","2","3","4","5")

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
hm1<-pheatmap(matrix, annotation_row = annot_row, annotation_colors = ann_colors, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,cluster_cols=F,cluster_rows=T,fontsize = 14, border_color = NA) 

ordered_matrix <- matrix[order(-matrix[,3],matrix[,5],-matrix[,8]),c(3,5,8)] #- Segerstolpe
#ordered_matrix <- matrix[order(-matrix[,3],matrix[,1],-matrix[,10],-matrix[,4]),c(3,1,10,4)] #- Baron
png("../data/Fig3d_2.png", units="in", width=7, height=7, res=300)
pheatmap(ordered_matrix,annotation_row = annot_row, annotation_colors = ann_colors, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,cluster_cols=F,cluster_rows=F,fontsize = 12, border_color = NA, cellwidth = 25)
dev.off()

annot_row$HbA1c = c(5.8,NA,5.4,NA,NA,5.5,NA,4.3,NA,NA,6.9,NA,5.4,5.8,5.5,4.6,6.1,5.3,4.5,NA,4.3,NA,5.4,
                    5.5,6.2,5.8,6.0,5.9,7.0,5.6,5.6,8.6,5.4,5.3,6.2,5.0,5.2,6.0,5.3,5.6,NA,6.8,6.8,5.5,5.3,8.0,
                    4.6,6.2,NA,5.6,6.0,5.5,6.2,6.4,7.0,10.0,5.1,5.3,5.4,5.3,6.0,5.8,5.4,5.7,6.2,5.2,5.7,5.2,5.4,
                    5.6,5.6,5.8,5.7,5.6,5.7,6.1,6.1,7.3,NA,5.7,5.5,7.2,6.4,6.7,5.9,5.6,6.0,5.6,5.3)

matf <- as.data.frame(mat)
matf$HbA1c <- annot_row$HbA1c
matf$Disease <- annot_row$Disease

meltedf<- melt(matf,id.vars = c("Disease", "HbA1c"))

meltedf <- na.omit(meltedf)
meltedf$Disease <- factor(meltedf$Disease, levels = c("N", "IGT", "T2D"))
png("../data/Fig3d_1.png", units="in", width=10, height=10, res=300)
ggscatter(meltedf, x="HbA1c", y="value",  color="Disease", palette = c("yellow","lightgreen","darkgreen"),
          facet.by = "variable", xlab = "HbA1c", ylab ="Cell type proportion", ncol = 5, xlim = c(5,10), ylim = c(0,0.8)) + stat_cor(method = "pearson", size = 3.2, label.x.npc = 0.2) +
          geom_smooth(color = "black", method=lm, se=TRUE, linetype="solid", fullrange=TRUE) + theme(title = element_text(size = 14, face = "bold"), text = element_text(size=12)) 
dev.off()

# Reordering Baron heatmap based on Segerstolpe heatmap
f <-read.csv(file="../data/Fig3_GEO_Fadista_5topics_metagene_normalized_Baron.csv",header=FALSE)
mat <- as.matrix(f[,c(1:9,12)])
baron_mat <- mat
colnames(mat)<-c("acinar","alpha","beta","delta","ductal","endothelial","epsilon","gamma","mast","PSC")

annot_row = data.frame(Disease=factor(c("N","na","N","na","na","N","na","N","na","na","T2D","na","N","N","N","N","IGT","N","N","na","N", "na", "N",
                                        "N","IGT","N", "IGT", "N",  "T2D", "N",  "N",  "T2D", "N", "N",  "IGT", "N", "N", "IGT",  "N",  "N",   "na",  "T2D",  "T2D",  "N",  "N",  "T2D",
                                        "N", "IGT",  "na", "N", "IGT",  "N", "IGT",  "IGT",  "T2D", "T2D",  "N",  "N",  "N",  "N",  "IGT",  "N",  "N",  "N",  "IGT",  "N",  "N",  "N",  "N",
                                        "N", "N",  "N",  "N",  "N",  "N",  "IGT",  "IGT",  "T2D",  "na",  "N",  "N",  "T2D",  "IGT",  "T2D",  "N",  "N",  "IGT",  "N",  "N")))

rownames(annot_row) <- rownames(f)
ann_colors = list(Disease=c(N = "yellow", IGT = "lightgreen", T2D = "darkgreen", na = "white" ),
                  HbA1c = c("orange","blue"))

matrix <-t(apply(mat,1, function(x) x/sum(x))) #Change each element by value / rowSum
rownames(matrix) <- rownames(f)
reordered_matrix <- matrix[hm1$tree_row$order,] #Re-order on the basis of segerstolpe
png("../data/Fig3e.png", units="in", width=7, height=7, res=300)
pheatmap(reordered_matrix, annotation_row = annot_row, annotation_colors = ann_colors, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,cluster_cols=F,cluster_rows=F,fontsize = 16, border_color = NA) 
dev.off()

cor.test(c(segerstolpe_mat),c(baron_mat))

annot_row$HbA1c = c(5.8,NA,5.4,NA,NA,5.5,NA,4.3,NA,NA,6.9,NA,5.4,5.8,5.5,4.6,6.1,5.3,4.5,NA,4.3,NA,5.4,
                    5.5,6.2,5.8,6.0,5.9,7.0,5.6,5.6,8.6,5.4,5.3,6.2,5.0,5.2,6.0,5.3,5.6,NA,6.8,6.8,5.5,5.3,8.0,
                    4.6,6.2,NA,5.6,6.0,5.5,6.2,6.4,7.0,10.0,5.1,5.3,5.4,5.3,6.0,5.8,5.4,5.7,6.2,5.2,5.7,5.2,5.4,
                    5.6,5.6,5.8,5.7,5.6,5.7,6.1,6.1,7.3,NA,5.7,5.5,7.2,6.4,6.7,5.9,5.6,6.0,5.6,5.3)

matf <- as.data.frame(mat)
matf$HbA1c <- annot_row$HbA1c
matf$Disease <- annot_row$Disease

meltedf<- melt(matf,id.vars = c("Disease", "HbA1c"))

meltedf <- na.omit(meltedf)
meltedf$Disease <- factor(meltedf$Disease, levels = c("N", "IGT", "T2D"))
png("../data/Fig3f_1.png", units="in", width=10, height=10, res=300)
ggscatter(meltedf, x="HbA1c", y="value",  color="Disease", palette = c("yellow","lightgreen","darkgreen"),
          facet.by = "variable", xlab = "HbA1c", ylab ="Cell type proportion", ncol = 2, xlim = c(5,10), ylim = c(0,0.8)) + stat_cor(method = "pearson", size = 3.2, label.x.npc = 0.2, label.y.npc = 1.0) +
geom_smooth(color = "black", method=lm, se=TRUE, linetype="solid", fullrange=TRUE) + theme(title = element_text(size = 16, face = "bold"), text = element_text(size=12)) 
dev.off()

ordered_matrix <- matrix[order(-matrix[,3],matrix[,1],-matrix[,8],-matrix[,4]),c(3,1,8,4)] #- Baron
png("../data/Fig3f_2.png", units="in", width=7, height=7, res=300)
pheatmap(ordered_matrix,annotation_row = annot_row, annotation_colors = ann_colors, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,cluster_cols=F,cluster_rows=F,fontsize = 12, border_color = NA, cellwidth = 25)
dev.off()