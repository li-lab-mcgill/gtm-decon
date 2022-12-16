library(pheatmap)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggpubr)


##Generate sparse matrix file from bulkRS files

f<-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/PAAD.pp.tab",sep="\t",row.names=1)
f<-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA.pp.tab",sep="\t",row.names=1)
#f <- f[intersect(names(f),coding_genes)]
#f2<-f[,-c(ncol(f)]
#f2<-f[,-c(20506,20505,20504,20503)] ## For BRCA.pp.TNBC_ER.tab - i dont know why
avg_val <- mean(apply(f,1,quantile,c(0.75)))
#avg_val <- mean(apply(f2,1,quantile,c(0.50)))
f[f<avg_val] <-0
#write.table(f,file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/PAAD.pp.sparse.75p.tab",sep="\t")
write.table(f,file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA.pp.sparse.75p.tab",sep="\t")

f <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/COMBINED/PANCREAS_IMMUNE/TCGA_PAAD_bulkData_trainData_JCVB0_nmar_K90_iter4_metagene_normalized.csv",header=FALSE)
mat <- as.matrix(f[,c(1:9,12,15:18)])

#mat <- as.matrix(f[,c(1:45,56:60)])
colnames(mat)<-c("acinar","alpha","beta","delta","ductal","endothelial","epsilon","gamma","mast","PSC","B_cell","Fibroblast","Monocyte","T_cell")
#colnames(mat)<-c("acinar","2","3","4","5","alpha","2","3","4","5","beta","2","3","4","5","delta","2","3","4","5","ductal","2","3","4","5","endothelial","2","3","4","5","epsilon","2","3","4","5","gamma","2","3","4","5","mast","2","3","4","5","PSC","2","3","4","5")


matrix <-t(apply(mat,1, function(x) x/sum(x))) #Change each element by value / rowSum
rownames(matrix) <- rownames(f)

pheatmap(matrix, annotation_row = annot_row, annotation_colors = ann_colors, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,cluster_cols=T,cluster_rows=T,fontsize = 12, border_color = NA) 
hm1<-pheatmap(matrix, annotation_row = annot_row, annotation_colors = ann_colors, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,cluster_cols=F,cluster_rows=T,fontsize = 14, border_color = NA) 

pheatmap(mat,color=colorRampPalette(c("white", "red"))(50),show_rownames = F,cluster_cols=F,cluster_rows=T,fontsize = 12, border_color = NA)


## Phi matrix ## 
f <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/BN/5topics_hvg/trainData_JCVB0_nmar_K35_iter4_phi_normalized.csv",header=FALSE)
genes <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/BN/5topics_hvg/genes.txt",header=FALSE)

f <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/BC_Wu_etal/pp/5topics/trainData_JCVB0_nmar_K145_iter4_phi_normalized.csv",header=FALSE)
genes <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/BC_Wu_etal/pp/5topics/genes.txt",header=FALSE)


df <- cbind(genes,f)
#rownames(df)<-df$V1

#Drop the first second and third column - they dont correspond to phi values
df<-df[-(2:3)]
#df<-df[-1]
#Get top50 genes per column
datalist = list()
topn_df <- data.frame()
for(i in 2:ncol(df)) {       # for-loop over columns
  subset_df <- df %>% arrange(desc(df[,i])) %>% slice(1:20)
  datalist[[i]]<-subset_df
}
topn_df <- do.call(rbind,datalist)
celltypes<-c("B cells Memory","B cells Naive","CAFs MSC iCAF-like","CAFs myCAF-like","Cancer Basal SC","Cancer Cycling","Cancer Her2 SC","Cancer LumA SC","Cancer LumB SC","Cycling_Myeloid","Cycling PVL","Cycling T-cells","DCs","Endothelial ACKR1","Endothelial CXCL12","Endothelial Lymphatic LYVE1","Endothelial RGS5","Luminal Progenitors","Macrophage","Mature Luminal","Monocyte","Myoepithelial","NK cells","NKT cells","Plasmablasts","PVL Differentiated","PVL Immature","T cells CD4+","T cells CD8+")
pheatmap(as.data.frame(select(topn_df,2:ncol(topn_df))),cluster_cols=F,cluster_rows=F,color=colorRampPalette(c("white","blue"))(50),show_rownames = F,fontsize = 10,labels_col = celltypes)
#pheatmap(as.data.frame(select(topn_df,2:ncol(topn_df))),cluster_cols=F,cluster_rows=F,color=colorRampPalette(c("white","blue"))(50),show_rownames = F,fontsize = 20,labels_col = celltypes)

# Get genesets for pancreatic tissue from CellMarker Database
cellmarker <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/Human_cell_markers_geneset.txt",sep="#")
cellMarker_gs <- split(x = cellmarker$GeneSymbol, f = cellmarker$TissueCellName)


# Pool all annotations for breast celltypes together from cellmarker db - as these are small genes
marker_list <- cellMarker_gs[["Breast - Basal epithelial cell"]]
a <- cellMarker_gs[["Breast - Basal progenitor cell"]]
marker_list <- append(marker_list,a)
a <- cellMarker_gs[["Breast - Cancer stem cell"]]
marker_list <- append(marker_list,a)
a <- cellMarker_gs[["Breast - Epithelial cell"]]
marker_list <- append(marker_list,a)
a <- cellMarker_gs[["Breast - Epithelial progenitor cell"]]
marker_list <- append(marker_list,a)
a <- cellMarker_gs[["Breast - Hematopoietic stem cell"]]
marker_list <- append(marker_list,a)
a <- cellMarker_gs[["Breast - Immune cell"]]
marker_list <- append(marker_list,a)
a <- cellMarker_gs[["Breast - Luminal epithelial cell"]]
marker_list <- append(marker_list,a)
a <- cellMarker_gs[["Breast - Luminal progenitor cell"]]
marker_list <- append(marker_list,a)
a <- cellMarker_gs[["Breast - Luminal progenitor"]]
marker_list <- append(marker_list,a)
a <- cellMarker_gs[["Breast - Mesenchymal stem cell"]]
marker_list <- append(marker_list,a)
a <- cellMarker_gs[["Breast - Myoepithelial cell"]]
marker_list <- append(marker_list,a)
a <- cellMarker_gs[["Breast - Natural killer cell"]]
marker_list <- append(marker_list,a)
a <- cellMarker_gs[["Breast - Progenitor cell"]]
marker_list <- append(marker_list,a)
cellMarker_gs[["Breast - all celltypes"]] <- marker_list

matrix <- as.data.frame(select(topn_df,2:ncol(topn_df)))
a <- topn_df$V1
annot_list_cellmarkerdb <- a %in% cellMarker_gs[["Breast - all celltypes"]]
sum(annot_list_cellmarkerdb,na.rm=TRUE)

mammary_epithelial_markers <- c("WNT5B","PRLR","CLDN4","CSN3","CSN1S1","RARRES1","KRT17","NNMT","ALDH1A3","SGMS2","PIEZO1","ANO1","KRT5","NEDD9","FOXP4","LSR","LIF","KRT19","KRT7","STAT3","KRT8","KRT18","BTN1A1","BTN2A2","BTN3A1","MMP14","TOX3","CDH1","FFAR1","PYGO2","SERPINB5","RUNX2","CSN2","CITED1","IRF6","SULT1E1","KRT14","KDM3A","SLC30A2","PADI2","CLDN1","ABCG2","PIP")
luminal_epithelial_markers <- c("AR","PIP","SERPINA1","ATP7B","HIF1A","FGFR4","DDR1","CEBPD","PGR","KRT8","UXT","PTH1R","CD9","AQP3","ATP2C2","WNT5A","SLC12A2","ESR1","AQP5","RUNX1","CDH1","MUC1","ANPEP","KLK3","KRT23","FGG","KRT18","ANKRD30A","SLPI","PROM1","KRT19","SYTL2","CD74","AGR2","LTF","SAA2","FGFR2","SERPINB4","SERPINB3","WFDC2","LCN2","BTG1","CLDN4","ANXA1","HMGA1","STC2","AREG","TNFSF10")

both_markers <- mammary_epithelial_markers
both_markers <- append(both_markers, luminal_epithelial_markers)
a <-topn_df$V1
annot_list_mammary <- a %in% mammary_epithelial_markers
a <-topn_df$V1
annot_list_luminal <- a %in% luminal_epithelial_markers
a <-topn_df$V1
annot_list_panglaodb <- a %in% both_markers

#annot_row = data.frame(Mammary_markers = factor(annot_list_mammary),
#                       Luminal_markers = factor(annot_list_luminal))
#ann_colors = list(Mammary_markers=c("TRUE" = "red", "FALSE" = "white", "NA" = "lightgrey"),
#                  Luminal_markers = c("TRUE" = "orange", "FALSE" = "white", "NA" = "lightgrey"))

annot_row = data.frame(CellMarkerDB = factor(annot_list_cellmarkerdb),
                       PanglaoDB = factor(annot_list_panglaodb))
ann_colors = list(CellMarkerDB=c("TRUE" = "red", "FALSE" = "white", "NA" = "lightgrey"),
                  PanglaoDB = c("TRUE" = "orange", "FALSE" = "white", "NA" = "lightgrey"))
rownames(topn_df) <- paste0("row_", seq(nrow(topn_df)))
rownames(annot_row) = rownames(topn_df)
matrix <- select(topn_df,2:26)
matrix <- matrix[1:500,]
max_val <- max(matrix)
matrix <- matrix/max_val
pheatmap(as.data.frame(matrix),annotation_row = annot_row, annotation_colors = ann_colors,cluster_cols=F,cluster_rows=F,color=colorRampPalette(c("white","blue"))(50),show_rownames = F,fontsize = 12,labels_col = c("Basal","Basal_Myoepithelial","Luminal_1_1","Luminal_1_2","Luminal_2"))


### GTEX - Breast tissue ###
f <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/COMBINED/BREAST_IMMUNE/GTEX_Breast_bulkData_combined_scaled_metagene_normalized.csv",header=FALSE)
f1 <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/COMBINED/BREAST_IMMUNE_hvg/GTEX_Breast_bulkData_combined_scaled_metagene_normalized.csv",header=FALSE)
f1 <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/BC_Wu_etal/pp/5topics/GTEX_Breast_bulkData_trainData_JCVB0_nmar_K145_iter4_metagene_normalized.csv",header=FALSE)
f1<-f1[-c(ncol(f1))]
#colnames(f1)=c("Basal","Basal_Myoepithelial","Luminal_1_1","Luminal_1_2","Luminal_2","Unclassified","Unnamed","B_cell","Fibroblast","Monocyte","T_cell")
colnames(f1)=c("B cells Memory","B cells Naive","CAFs MSC iCAF-like","CAFs myCAF-like","Cancer Basal SC","Cancer Cycling","Cancer Her2 SC","Cancer LumA SC","Cancer LumB SC","Cycling_Myeloid","Cycling PVL","Cycling T-cells","DCs","Endothelial ACKR1","Endothelial CXCL12","Endothelial Lymphatic LYVE1","Endothelial RGS5","Luminal Progenitors","Macrophage","Mature Luminal","Monocyte","Myoepithelial","NK cells","NKT cells","Plasmablasts","PVL Differentiated","PVL Immature","T cells CD4+","T cells CD8+")
#annot_col = data.frame(cell_type = factor(rep(c("Breast","Immune"), c(,4))))
#rownames(annot_col) = colnames(f1)
pheatmap(f1,show_rownames = F,color=colorRampPalette(c("white", "red"))(50),fontsize = 12, fontsize_col = 12,cellwidth=20,cluster_cols = F)
pheatmap(f1,show_rownames = F,annotation = annot_col,color=colorRampPalette(c("white", "red"))(50),fontsize = 12, fontsize_col = 12,cellwidth=20,cluster_cols = F)

## TCGA - BRCA ##
f2 <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/COMBINED/BREAST_IMMUNE_hvg/TCGA_BRCA_bulkData_combined_scaled_metagene_normalized.csv",header=FALSE)
f2 <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/BC_Wu_etal/pp/5topics/TCGA_BRCA_bulkData_trainData_JCVB0_nmar_K145_iter4_metagene_normalized.csv",header=FALSE)

f2<-f2[-c(ncol(f2))]
colnames(f2)=c("B cells Memory","B cells Naive","CAFs MSC iCAF-like","CAFs myCAF-like","Cancer Basal SC","Cancer Cycling","Cancer Her2 SC","Cancer LumA SC","Cancer LumB SC","Cycling_Myeloid","Cycling PVL","Cycling T-cells","DCs","Endothelial ACKR1","Endothelial CXCL12","Endothelial Lymphatic LYVE1","Endothelial RGS5","Luminal Progenitors","Macrophage","Mature Luminal","Monocyte","Myoepithelial","NK cells","NKT cells","Plasmablasts","PVL Differentiated","PVL Immature","T cells CD4+","T cells CD8+")
#colnames(f2)=c("Basal","Basal_Myoepithelial","Luminal_1_1","Luminal_1_2","Luminal_2","Unclassified","Unnamed","B_cell","Fibroblast","Monocyte","T_cell")
#annot_col = data.frame(cell_type = factor(rep(c("Breast","Immune"), c(7,4))))
#rownames(annot_col) = colnames(f2)
pheatmap(f2,show_rownames = F,color=colorRampPalette(c("white", "red"))(50),fontsize = 12, fontsize_col = 12,cellwidth=20,cluster_cols = F)
pheatmap(f2,show_rownames = F,annotation = annot_col,color=colorRampPalette(c("white", "red"))(50),fontsize = 12, fontsize_col = 12,cellwidth=20,cluster_cols = F)

gtex_tcga <- rbind(f1,f2)
rownames(gtex_tcga) <- paste0("row_", seq(nrow(gtex_tcga)))

pheatmap(gtex_tcga)
labels <- rep(c("GTEx","TCGA"), c(nrow(f1),nrow(f2)))
annot_row <- data.frame(DataSource = labels)
row.names(annot_row) <- row.names(gtex_tcga)
annot_colors = list(DataSource = c(GTEx="purple",TCGA="darkgreen"))
pheatmap(gtex_tcga, annotation_row = annot_row, annotation_colors = annot_colors,color=colorRampPalette(c("white", "red"))(50),show_rownames=F,fontsize=12) 
pheatmap(gtex_tcga[c(1:5,8:11)], annotation_row = annot_row, annotation_colors = annot_colors) 
pheatmap(gtex_tcga[c(1:5,8:11)], show_rownames = F,annotation_row = annot_row, annotation_colors = annot_colors, color=colorRampPalette(c("white", "red"))(50),fontsize = 12, cellwidth=20)
pheatmap(gtex_tcga[c(1:5)], show_rownames = F,annotation_row = annot_row, annotation_colors = annot_colors, color=colorRampPalette(c("white", "red"))(50),fontsize = 12, cellwidth=20)

## Heatmap - celltype proportions from Wu etal. scRNAseq
f1 <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts_gtm-decon/BC_Wu_etal/celltype_prop_minor_bRS_norm.csv",row.names=1)
#f2<-f2[-c(ncol(f2))]
#colnames(f2)=c("B cells Memory","B cells Naive","CAFs MSC iCAF-like","CAFs myCAF-like","Cancer Basal SC","Cancer Cycling","Cancer Her2 SC","Cancer LumA SC","Cancer LumB SC","Cycling_Myeloid","Cycling PVL","Cycling T-cells","DCs","Endothelial ACKR1","Endothelial CXCL12","Endothelial Lymphatic LYVE1","Endothelial RGS5","Luminal Progenitors","Macrophage","Mature Luminal","Monocyte","Myoepithelial","NK cells","NKT cells","Plasmablasts","PVL Differentiated","PVL Immature","T cells CD4+","T cells CD8+")
#colnames(f2)=c("Basal","Basal_Myoepithelial","Luminal_1_1","Luminal_1_2","Luminal_2","Unclassified","Unnamed","B_cell","Fibroblast","Monocyte","T_cell")
#annot_col = data.frame(cell_type = factor(rep(c("Breast","Immune"), c(7,4))))
#rownames(annot_col) = colnames(f2)
pheatmap(f1,show_rownames = F,color=colorRampPalette(c("white", "red"))(50),fontsize = 12, fontsize_col = 12,cellwidth=20,cluster_cols = F,cluster_rows=F)

rownames(f1)<-c()
colnames(f1)<-c()

## Heatmap - bRS from Wu_etal
f2 <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/BC_Wu_etal/pp/5topics/breast_cancer_Wu_etal_Breast_bulkData_trainData_JCVB0_nmar_K145_iter4_metagene_normalized.csv",header=FALSE)
f2<-f2[-c(ncol(f2))]
colnames(f2)=c("B cells Memory","B cells Naive","CAFs MSC iCAF-like","CAFs myCAF-like","Cancer Basal SC","Cancer Cycling","Cancer Her2 SC","Cancer LumA SC","Cancer LumB SC","Cycling_Myeloid","Cycling PVL","Cycling T-cells","DCs","Endothelial ACKR1","Endothelial CXCL12","Endothelial Lymphatic LYVE1","Endothelial RGS5","Luminal Progenitors","Macrophage","Mature Luminal","Monocyte","Myoepithelial","NK cells","NKT cells","Plasmablasts","PVL Differentiated","PVL Immature","T cells CD4+","T cells CD8+")
#colnames(f2)=c("Basal","Basal_Myoepithelial","Luminal_1_1","Luminal_1_2","Luminal_2","Unclassified","Unnamed","B_cell","Fibroblast","Monocyte","T_cell")
#annot_col = data.frame(cell_type = factor(rep(c("Breast","Immune"), c(7,4))))
#rownames(annot_col) = colnames(f2)
pheatmap(f2,show_rownames = F,color=colorRampPalette(c("white", "red"))(50),fontsize = 12, fontsize_col = 12,cellwidth=20,cluster_cols = F,cluster_rows=F)

rownames(f2)<-c()
colnames(f2)<-c()


## Heatmap
## BRCA - extra ##
clinical <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts_gtm-decon/COMBINED/BRCA_clinical_notes_extra_breast_immune_scaled_10topics_UNIQ.csv")
clinical <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts_gtm-decon/COMBINED/BRCA_clinical_notes_extra_breast_immune_scaled_10topics_hvg_UNIQ.csv")
#clinical <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts_gtm-decon/COMBINED/PAAD_clinical_notes_pancreas_immune_combined_unscaled_i4_k10_UNIQ.csv")
clinical <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts_gtm-decon/BC_Wu_etal/BRCA_clinical_notes_extra_BC_Wu_etal_iter4_10topics_pp_UNIQ.csv")
annot_row <- data.frame("subtype"=clinical$histological_type)
row.names(annot_row) <- row.names(clinical)
#annot_row$gender <- clinical$gender
annot_row$stage <- clinical$pathologic_stage
annot_row$days_to_death <- as.integer(clinical$days_to_death)
annot_row$PAM50 <- clinical$PAM50
annot_row$ER_IHC <- clinical$ER_IHC
#annot_row$PR_IHC <- clinical$PR_IHC
#annot_row$HER2_IHC <- clinical$HER2_IHC
#annot_row$TumorPurity <- as.double(clinical$TumorPurity)
annot_row$ProliferationScore <- as.double(clinical$ProliferationScore)

#annot_col = data.frame(cell_type = factor(rep(c("Clinical","Breast","Immune","Topic"), c(30,7,4,10))))
#annot_col = data.frame(cell_type = factor(rep(c("Clinical","Breast","Immune"), c(18,7,25))))
annot_col = data.frame(cell_type = factor(rep(c("Clinical","Breast","Topic"), c(30,29,10))))
rownames(annot_col) = colnames(clinical)

#Set2,,Set2,,,YlGnBu,
annot_colors = list(
  subtype=c(infiltrating_lobular_carcinoma = "#66C2A5", infiltrating_ductal_carcinoma = "#808080", infiltrating_carcinoma_nos = "#8DA0CB", mucinous_carcinoma = "#E78AC3", medullary_carcinoma = "#A6D854", metaplastic_carcinoma = "#FFD92F", mixed_histology = "#E5C494", other = "#B3B3B3", na = "#FFFFFF"),
  gender=c(female = "pink", male = "blue"),
  stage=c(stage_i = "#66C2A5",stage_ii = "#FC8D62",stage_iii = "#8DA0CB",stage_iv = "#E78AC3",stage_x = "#A6D854", na = "#FFFFFF"),
  T_stage=c(t1 = "#FFFFCC",t2 = "#FED976",t3 = "#FD8D3C",t4 = "#E31A1C",tx = "#FFFFFF", na="#FFFFFF"),
  M_stage=c(m0 = "#FFFFFF", m1 ="#969696",mx = "#FFFFFF"),
  N_stage=c(n0 = "#FFFFD9",n1 = "#C7E9B4", n2 = "#41B6C4", n3 ="#225EA8", nx = "#FFFFFF", na="#FFFFFF"),
  days_to_death=c("darkgreen","white"),
  PAM50=c(Basal = "red", Her2 = "green", LumA = "lightblue", LumB = "darkblue" , Normal = "yellow", na = "white"),
  ER_IHC=c(Positive = "darkgreen", Negative = "red", Indeterminate = "grey", na = "white", Not_Evaluated = "lightgrey"),
  PR_IHC=c(Positive = "darkgreen", Negative = "red", Indeterminate = "grey", na = "white", Not_Evaluated = "lightgrey"),
  HER2_IHC=c(Positive = "darkgreen", Negative = "red", Indeterminate = "grey", na = "white", Not_Evaluated = "lightgrey", Not_Available = "white", Equivocal = "yellow"),
  TumorPurity=c("white","darkgrey"),
  ProliferationScore=c("red","white","blue"),
  cell_type=c(Breast = "#0000FF", Immune = "#00FF00", Topic = "#A9A9A9")
)

#Set1 - all topics retained
#m <- as.data.frame(select(clinical,31:35,38:ncol(clinical)))
m <- as.data.frame(select(clinical,31:ncol(clinical)))
#m <- as.data.frame(select(clinical,19:23))
matrix <-t(apply(m,1, function(x) x/sum(x))) #Change each element by value / rowSum
rownames(matrix)<-rownames(clinical)
pheatmap(matrix[,1:29],annotation_row = annot_row, annotation = annot_col, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 8, fontsize_col = 12,cellwidth = 20)
pheatmap(matrix,annotation_row = annot_row, annotation = annot_col, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 8, fontsize_col = 12,cellwidth = 20)


## Order based on fibroblast proportion
ordered_matrix<-matrix[order(-matrix[,5]),]
pheatmap(ordered_matrix,annotation_row = annot_row, annotation = annot_col, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 10,cluster_rows=F,cellwidth = 20)

## Correlation with stage or other params
subset_matrix <- clinical[,c('Basal','Basal_Myoepithelial','Luminal_1_1','Luminal_1_2','Luminal_2','Fibroblast')]
subset_matrix <- clinical[,c('B.cells.Memory','B.cells.Naive','CAFs.MSC.iCAF.like','CAFs.myCAF.like','Cancer.Basal.SC','Cancer.Cycling','Cancer.Her2.SC','Cancer.LumA.SC','Cancer.LumB.SC','Cycling_Myeloid','Cycling.PVL','Cycling.T.cells','DCs','Endothelial.ACKR1','Endothelial.CXCL12','Endothelial.Lymphatic.LYVE1','Endothelial.RGS5','Luminal.Progenitors','Macrophage','Mature.Luminal','Monocyte','Myoepithelial','NK.cells','NKT.cells','Plasmablasts','PVL.Differentiated','PVL.Immature','T.cells.CD4.','T.cells.CD8.')]
matf <- as.data.frame(subset_matrix)
matf$days_to_death <- clinical$days_to_death
matf$stage <- clinical$pathologic_stage

matf$stage <- gsub("ia","i",matf$stage)
matf$stage <- gsub("ib","i",matf$stage)

matf$stage <- factor(matf$stage,levels = c("stage_i","stage_ii","stage_iii","stage_iv","stage_X","na"))

matf$T_stage <- clinical$pathology_T_stage
matf$M_stage <- clinical$pathology_M_stage
matf$N_stage <- clinical$pathology_N_stage
matf$T_stage <- factor(matf$T_stage,levels = c("t1","t1b","t1c","t2","t3","t4","t4b"))
matf$N_stage <- factor(matf$N_stage,levels = c("n0","n0 (i+)","n0 (i-)","n0 (mol+)","n1","n1a","n1b","n1mi","n2","n2a","n3","n3a","n3b","nx"))                   
matf$M_stage <- factor(matf$M_stage,levels = c("m0","cm0 (i+)","m1","mx"))

matf$lymph_nodes <- clinical$number_of_lymph_nodes
matf$PAM50 <- clinical$PAM50
matf$PAM50 <- factor(matf$PAM50, levels = c("Normal","Basal","LumA","LumB","Her2","na"))
matf$ProliferationScore <- clinical$ProliferationScore

meltedf<- melt(matf,id.vars = c("stage", "days_to_death","T_stage","N_stage","M_stage","lymph_nodes","PAM50","ProliferationScore"))
meltedf <- na.omit(meltedf)

## Comparison with PAM50
subsetf <- meltedf[meltedf$PAM50 != "na",]

compare_means(value ~ PAM50, data = smallset)
my_comparisons <- list(c("Normal","Basal"),c("Normal","LumA"), c("Normal","LumB"), c("Normal","Her2"), c("Basal", "LumA"), c("Basal", "LumB"), c("Basal", "Her2"), c("LumA","LumB"), c("LumA","Her2"), c("LumB","Her2"))
ggboxplot(subsetf,x="PAM50", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "PAM50 classification", ylab = "Cell-type proportion") + font("xy.text",size = 14) + font ("ylab", size = 14) +
    stat_compare_means()
    stat_compare_means(comparisons = my_comparisons)
subsetf <- subsetf[subsetf$variable == "Basal" | subsetf$variable == "Luminal_1_1" | subsetf$variable == "Luminal_2",]

subsetf <- subsetf[subsetf$variable == "Cancer.Basal.SC" | subsetf$variable == "Cancer.Cycling" | subsetf$variable == "Cancer.Her2.SC" | subsetf$variable == "Cancer.LumA.SC" | subsetf$variable == "Cancer.LumB.SC" |subsetf$variable == "Endothelial.ACKR1" | subsetf$variable == "Luminal.Progenitors" | subsetf$variable == "Myoepithelial",]
ggboxplot(subsetf, x="PAM50", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "PAM50 classification", ylab = "Cell-type proportion", ncol=2) + font("xy.text",size = 10) + font ("ylab", size = 12) +
  stat_compare_means(position = position_nudge(x=0.5,y=-0.05))
  stat_compare_means(method="t.test") 
ggboxplot(meltedf, x="PAM50", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "PAM50", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

## Comparison with ProliferationScore
ggboxplot(meltedf, x="ProliferationScore", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "ProliferationScore", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 
ggscatter(meltedf, x="ProliferationScore", y="value",  color="PAM50", facet.by = "variable") + stat_compare_means()

## Comparison with stage
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "stage", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

ggboxplot(meltedf, x="stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

## Comparison with T_stage
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "T_stage", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

subsetf <- subsetf[subsetf$variable == "CAFs.MSC.iCAF.like" | subsetf$variable == "CAFs.myCAF.like" | subsetf$variable == "Cancer.cycling" | subsetf$variable == "Cancer.LumA.SC" | subsetf$variable == "Cancer.LumB.SC",]
ggboxplot(subsetf, x="T_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 8) +
  stat_compare_means() 

subsetf <- meltedf[meltedf$variable == "Ductal",]
ggboxplot(subsetf, x="T_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 14) + font("ylab",size=14) +
  stat_compare_means() 


## Comparison with N_stage
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "N_stage", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

ggboxplot(meltedf, x="N_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 


## Comparison with M_stage

ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "M_stage", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means(comparisons = my_comparisons) 

ggboxplot(meltedf, x="M_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

subsetf <- meltedf[meltedf$variable == "Ductal",]
ggboxplot(subsetf, x="M_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 14) + font("ylab",size=14) +
  stat_compare_means() 

subsetf <- meltedf[meltedf$variable == "Alpha",]
ggboxplot(subsetf, x="M_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

## COmparison with number of lymph nodes
ggboxplot(meltedf, x="lymph_nodes", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

ggboxplot(meltedf, x="days_to_death", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

## Correlation with stage or other params - for topics
subset_matrix <- clinical[,c('Topic1','Topic2','Topic3','Topic4','Topic5','Topic6','Topic7','Topic8','Topic9','Topic10')]
matf <- as.data.frame(subset_matrix)
matf$days_to_death <- clinical$days_to_death
matf$stage <- clinical$pathologic_stage

matf$stage <- gsub("ia","i",matf$stage)
matf$stage <- gsub("ib","i",matf$stage)

matf$stage <- factor(matf$stage,levels = c("stage_i","stage_ii","stage_iii","stage_iv","stage_X","na"))

matf$T_stage <- clinical$pathology_T_stage
matf$M_stage <- clinical$pathology_M_stage
matf$N_stage <- clinical$pathology_N_stage
matf$T_stage <- factor(matf$T_stage,levels = c("t1","t1b","t1c","t2","t3","t4","t4b"))
matf$N_stage <- factor(matf$N_stage,levels = c("n0","n0 (i+)","n0 (i-)","n0 (mol+)","n1","n1a","n1b","n1mi","n2","n2a","n3","n3a","n3b","nx"))                   
matf$M_stage <- factor(matf$M_stage,levels = c("m0","cm0 (i+)","m1","mx"))

matf$lymph_nodes <- clinical$number_of_lymph_nodes
matf$PAM50 <- clinical$PAM50
matf$PAM50 <- factor(matf$PAM50, levels = c("Normal","Basal","LumA","LumB","Her2","na"))
matf$ProliferationScore <- clinical$ProliferationScore

meltedf<- melt(matf,id.vars = c("stage", "days_to_death","T_stage","N_stage","M_stage","lymph_nodes","PAM50","ProliferationScore"))
meltedf <- na.omit(meltedf)

## Comparison with PAM50
subsetf <- meltedf[meltedf$PAM50 != "na",]


subsetf <- subsetf[subsetf$variable == "Topic7" | subsetf$variable == "Topic8" | subsetf$variable == "Topic5",]
my_comparisons <- list(c("Normal","Basal"),c("Normal","LumA"), c("Normal","LumB"), c("Normal","Her2"), c("Basal", "LumA"), c("Basal", "LumB"), c("Basal", "Her2"), c("LumA","LumB"), c("LumA","Her2"), c("LumB","Her2"))
ggboxplot(subsetf,x="PAM50", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "PAM50 classification", ylab = "Cell-type proportion", ncol=4) + font("xy.text",size = 14) + font ("ylab", size = 14) +
  stat_compare_means(comparisons = my_comparisons)

ggboxplot(subsetf, x="PAM50", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "PAM50 classification", ylab = "Cell-type proportion", ncol=4) + font("xy.text",size = 14) + font("ylab",size = 14) +
  stat_compare_means() 

ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "PAM50", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

## Comparison with ProliferationScore
ggboxplot(meltedf, x="ProliferationScore", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "ProliferationScore", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 
ggscatter(meltedf, x="ProliferationScore", y="value",  color="PAM50", facet.by = "variable")

## Comparison with stage
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "stage", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

ggboxplot(meltedf, x="stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

## Comparison with T_stage
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "T_stage", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

ggboxplot(meltedf, x="T_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

subsetf <- meltedf[meltedf$variable == "Ductal",]
ggboxplot(subsetf, x="T_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 14) + font("ylab",size=14) +
  stat_compare_means() 


## Comparison with N_stage
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "N_stage", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

ggboxplot(meltedf, x="N_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 


## Comparison with M_stage
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "M_stage", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

ggboxplot(meltedf, x="M_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

subsetf <- meltedf[meltedf$variable == "Ductal",]
ggboxplot(subsetf, x="M_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 14) + font("ylab",size=14) +
  stat_compare_means() 

subsetf <- meltedf[meltedf$variable == "Alpha",]
ggboxplot(subsetf, x="M_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

## COmparison with number of lymph nodes
ggboxplot(meltedf, x="lymph_nodes", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

ggboxplot(meltedf, x="days_to_death", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 



##Tissue proportion vs. Immune proportion
ti <- as.data.frame(select(clinical,31:35))
im <- as.data.frame(select(clinical,38:41))

tissue_prop <- rowSums(ti)
immune_prop <- rowSums(im)

tissue_immune <- data.frame(tissue_prop,immune_prop)
rownames(tissue_immune) <- rownames(clinical)
rownames(annot_row) <- rownames(clinical)


pheatmap(tissue_immune[order(tissue_immune[,2]),],annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 10,cluster_rows = F)

pheatmap(tissue_immune,annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 10)

annot_row$lymph_nodes <- clinical$number_of_lymph_nodes
annot_row$T_stage <- clinical$pathology_T_stage
annot_row$M_stage <- clinical$pathology_M_stage
annot_row$N_stage <- clinical$pathology_N_stage

annot_row$T_stage <- gsub("^$","na",annot_row$T_stage)
annot_row$M_stage <- gsub("^$","na",annot_row$M_stage)
annot_row$N_stage <- gsub("^$","na",annot_row$N_stage)


matf <- as.data.frame(tissue_immune)
matf$days_to_death <- annot_row$days_to_death
matf$stage <- annot_row$stage
matf$PAM50 <- annot_row$PAM50
matf$subtype <- annot_row$subtype
#matf$fibroblast <- matrix[,c(7)]

matf$stage <- gsub("ia","i",matf$stage)
matf$stage <- gsub("ib","i",matf$stage)

matf$T_stage <- annot_row$T_stage
matf$M_stage <- annot_row$M_stage
matf$N_stage <- annot_row$N_stage
matf$lymph_nodes <- annot_row$lymph_nodes

meltedf<- melt(matf,id.vars = c("stage", "days_to_death","T_stage","N_stage","M_stage","lymph_nodes","PAM50","subtype"))
meltedf <- na.omit(meltedf)
ggscatter(meltedf, x="days_to_death", y="value",  color="stage", facet.by = "variable")
ggscatter(meltedf, x="days_to_death", y="value",  color="M_stage", facet.by = "variable")
ggscatter(meltedf, x="days_to_death", y="value",  color="T_stage", facet.by = "variable")
ggscatter(meltedf, x="days_to_death", y="value",  color="N_stage", facet.by = "variable")

ggscatter(meltedf, x="lymph_nodes", y="value",  color="stage", facet.by = "variable")
ggscatter(meltedf, x="lymph_nodes", y="value",  color="M_stage", facet.by = "variable")
ggscatter(meltedf, x="lymph_nodes", y="value",  color="T_stage", facet.by = "variable")
ggscatter(meltedf, x="lymph_nodes", y="value",  color="N_stage", facet.by = "variable")

meltedf <- meltedf[meltedf$stage!="na",]
meltedf <- meltedf[meltedf$stage!="stage_x",]
ggboxplot(meltedf, x="variable", y="value",  color="variable", palette =c("#00AFBB", "#FC4E07"), add = "jitter", facet.by = "stage", xlab = "", ylab = "Combined cell-type proportions per tissue", ncol=4) +
  stat_compare_means()

ggboxplot(meltedf, x="variable", y="value",  color="variable", palette =c("#00AFBB", "#FC4E07"), add = "jitter", facet.by = "T_stage", xlab = "", ylab = "Combined cell-type proportions per tissue") +
  stat_compare_means()

ggboxplot(meltedf, x="variable", y="value",  color="variable", palette =c("#00AFBB", "#FC4E07"), add = "jitter", facet.by = "M_stage", xlab = "", ylab = "Combined cell-type proportions per tissue") +
  stat_compare_means()

ggboxplot(meltedf, x="variable", y="value",  color="variable", palette =c("#00AFBB", "#FC4E07"), add = "jitter", facet.by = "N_stage", xlab = "", ylab = "Combined cell-type proportions per tissue") +
  stat_compare_means()

meltedf <- meltedf[meltedf$PAM50!="na",]
ggboxplot(meltedf, x="variable", y="value",  color="variable", palette =c("#00AFBB", "#FC4E07"), add = "jitter", facet.by = "PAM50", xlab = "", ylab = "Combined cell-type proportions per tissue", ncol=5) +
  stat_compare_means()

ggboxplot(meltedf, x="variable", y="value",  color="variable", palette =c("#00AFBB", "#FC4E07"), add = "jitter", facet.by = "subtype", xlab = "", ylab = "Combined cell-type proportions per tissue") +
  stat_compare_means()


meltedf_f <- meltedf[meltedf$variable=='fibroblast',]
ggboxplot(meltedf_f, x="stage", y="value",  color="stage", palette =c("yellow","orange","red","brown"), add = "jitter",  xlab = "", ylab = "Proportion of fibroblast") +
  stat_compare_means()


# b) Kaplan - Meier curves for TCGA PAAD dataset

library("survival")
library("survminer")
library("cluster")
library("factoextra")

#Breast cancer - with immune dataset
#d <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts_gtm-decon/COMBINED/BRCA_clinical_notes_extra_breast_immune_scaled_10topics_hvg_UNIQ.csv")
d <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts_gtm-decon/BC_Wu_etal/BRCA_clinical_notes_extra_BC_Wu_etal_iter4_10topics_pp_UNIQ.csv")

d$days_to_death <- as.numeric(d$days_to_death)
d$days_to_last_followup <- as.numeric(d$days_to_last_followup)
#d <- d[d$days_to_death != "na",]

#'B.cells.Memory','B.cells.Naive','CAFs.MSC.iCAF.like','CAFs.myCAF.like','Cancer.Basal.SC','Cancer.Cycling','Cancer.Her2.SC','Cancer.LumA.SC','Cancer.LumB.SC','Cycling_Myeloid','Cycling.PVL','Cycling.T.cells','DCs','Endothelial.ACKR1','Endothelial.CXCL12','Endothelial.Lymphatic.LYVE1','Endothelial.RGS5','Luminal.Progenitors','Macrophage','Mature.Luminal','Monocyte','Myoepithelial','NK.cells','NKT.cells','Plasmablasts','PVL.Differentiated','PVL.Immature','T.cells.CD4.','T.cells.CD8.')
covariates <- c("years_to_birth", "gender",  "pathologic_stage", "pathology_T_stage", "pathology_N_stage",
                "pathology_M_stage","histological_type","number_of_lymph_nodes",
                "PAM50","ER_IHC","PR_IHC","HER2_IHC","ProliferationScore",
                "B.cells.Memory","B.cells.Naive","CAFs.MSC.iCAF.like","CAFs.myCAF.like","Cancer.Basal.SC","Cancer.Cycling","Cancer.Her2.SC","Cancer.LumA.SC","Cancer.LumB.SC","Cycling_Myeloid","Cycling.PVL","Cycling.T.cells","DCs","Endothelial.ACKR1","Endothelial.CXCL12","Endothelial.Lymphatic.LYVE1","Endothelial.RGS5","Luminal.Progenitors","Macrophage","Mature.Luminal","Monocyte","Myoepithelial","NK.cells","NKT.cells","Plasmablasts","PVL.Differentiated","PVL.Immature","T.cells.CD4.","T.cells.CD8.",
                "Topic1","Topic2","Topic3","Topic4","Topic5","Topic6","Topic7",
                "Topic8","Topic9","Topic10")

covariates <- c("Basal","Basal_Myoepithelial","Luminal_1_1","Luminal_1_2","Luminal_2","Unclassified","Unnamed",
                "B_cell","Fibroblast","Monocyte","T_cell","Topic1","Topic2","Topic3","Topic4","Topic5","Topic6","Topic7",
                "Topic8","Topic9","Topic10")

d$time <- ifelse(d$vital_status == 1,d$days_to_death,d$days_to_last_followup)

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, vital_status)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = d)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         Delta<-signif(x$coef[1], digits=2);#coeficient Delta
                         HR <-signif(x$coef[2], digits=2);#exp(Delta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(Delta, HR, wald.test, p.value)
                         names(res)<-c("Delta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

# Multivariate cox regression
res.cox <- coxph(Surv(time, vital_status) ~ Basal + Basal_Myoepithelial + Luminal_1_1 + Luminal_1_2 + Luminal_2 + Unclassified +
                   B_cell + Fibroblast + Monocyte + T_cell + Topic1 + Topic2 + Topic3 + Topic4 + Topic5 + Topic6 + Topic7 +
                   Topic8 + Topic9 + Topic10 , data = d)

res.cox <- coxph(Surv(time+ vital_status) ~ B.cells.Memory + B.cells.Naive + CAFs.MSC.iCAF.like + CAFs.myCAF.like + Cancer.Basal.SC+Cancer.Cycling + Cancer.Her2.SC+Cancer.LumA.SC + Cancer.LumB.SC+Cycling_Myeloid + Cycling.PVL+Cycling.T.cells + DCs + Endothelial.ACKR1 + Endothelial.CXCL12 + Endothelial.Lymphatic.LYVE1 + Endothelial.RGS5 + Luminal.Progenitors + Macrophage + Mature.Luminal + Monocyte + Myoepithelial + NK.cells + NKT.cells + Plasmablasts + PVL.Differentiated + PVL.Immature + T.cells.CD4. + T.cells.CD8. +
                 Topic1 + Topic2 + Topic3 + Topic4 + Topic5 + Topic6 + Topic7 +
                 Topic8 + Topic9 + Topic10 , data = d)

summary(res.cox)


# Survival curves
#Access to values returned by survfit
#v <- data.frame(time = fit$time,
#                n.risk = fit$n.risk,
#                n.event = fit$n.event,
#                n.censor = fit$n.censor,
#                surv = fit$surv,
#                upper = fit$upper,
#                lower = fit$lower
#)
#head(v)

##################### Survival analysis on proportion of B_cell ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$B_cell), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$B_cell), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$B_cell, 2)
# get cluster means
aggregate(d$B_cell,by=list(prop_cluster$cluster),FUN=mean)
d$B_cell_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ B_cell_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$B_cell_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of T_cell ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$T_cell), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$T_cell), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$T_cell, 2)
# get cluster means
aggregate(d$T_cell,by=list(prop_cluster$cluster),FUN=mean)
d$T_cell_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ T_cell_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$T_cell_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Monocyte ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Monocyte), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Monocyte), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Monocyte, 2)
# get cluster means
aggregate(d$Monocyte,by=list(prop_cluster$cluster),FUN=mean)
d$Monocyte_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Monocyte_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Monocyte_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Fibroblast ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Fibroblast), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Fibroblast), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Fibroblast, 2)
# get cluster means
aggregate(d$Fibroblast,by=list(prop_cluster$cluster),FUN=mean)
d$Fibroblast_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Fibroblast_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Fibroblast_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


##################### Survival analysis on proportion of Basal ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Basal), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Basal), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Basal, 2)
# get cluster means
aggregate(d$Basal,by=list(prop_cluster$cluster),FUN=mean)
d$Basal_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Basal_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Basal_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


##################### Survival analysis on proportion of Basal_Myoepithelial ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Basal_Myoepithelial), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Basal_Myoepithelial), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Basal_Myoepithelial, 2)
# get cluster means
aggregate(d$Basal_Myoepithelial,by=list(prop_cluster$cluster),FUN=mean)
d$Basal_Myoepithelial_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Basal_Myoepithelial_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Basal_Myoepithelial_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


##################### Survival analysis on proportion of Luminal_1_1 ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Luminal_1_1), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Luminal_1_1), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Luminal_1_1, 2)
# get cluster means
aggregate(d$Luminal_1_1,by=list(prop_cluster$cluster),FUN=mean)
d$Luminal_1_1_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Luminal_1_1_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Luminal_1_1_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Luminal_1_2 ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Luminal_1_2), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Luminal_1_2), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Luminal_1_2, 2)
# get cluster means
aggregate(d$Luminal_1_2,by=list(prop_cluster$cluster),FUN=mean)
d$Luminal_1_2_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Luminal_1_2_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Luminal_1_2_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Luminal_2 ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Luminal_2), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Luminal_2), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Luminal_2, 2)
# get cluster means
aggregate(d$Luminal_2,by=list(prop_cluster$cluster),FUN=mean)
d$Luminal_2_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Luminal_2_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Luminal_2_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


##################### Survival analysis on proportion of Topic ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Topic9), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Topic9), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Topic9, 2)
# get cluster means
aggregate(d$Topic9,by=list(prop_cluster$cluster),FUN=mean)
d$Topic9_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Topic9_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Topic9_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of B.cells.Memory ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$B.cells.Memory), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$B.cells.Memory), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$B.cells.Memory, 2)
# get cluster means
aggregate(d$B.cells.Memory,by=list(prop_cluster$cluster),FUN=mean)
d$B.cells.Memory_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ B.cells.Memory_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$B.cells.Memory_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of B.cells.Memory ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$B.cells.Memory), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$B.cells.Memory), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$B.cells.Memory, 2)
# get cluster means
aggregate(d$B.cells.Memory,by=list(prop_cluster$cluster),FUN=mean)
d$B.cells.Memory_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ B.cells.Memory_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$B.cells.Memory_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


##################### Survival analysis on proportion of B.cells.Naive ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$B.cells.Naive), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$B.cells.Naive), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$B.cells.Naive, 2)
# get cluster means
aggregate(d$B.cells.Naive,by=list(prop_cluster$cluster),FUN=mean)
d$B.cells.Naive_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ B.cells.Naive_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$B.cells.Naive_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of CAFs.MSC.iCAF.like ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$CAFs.MSC.iCAF.like), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$CAFs.MSC.iCAF.like), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$CAFs.MSC.iCAF.like, 2)
# get cluster means
aggregate(d$CAFs.MSC.iCAF.like,by=list(prop_cluster$cluster),FUN=mean)
d$CAFs.MSC.iCAF.like_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ CAFs.MSC.iCAF.like_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$CAFs.MSC.iCAF.like_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of CAFs.myCAF.like ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$CAFs.myCAF.like), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$CAFs.myCAF.like), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$CAFs.myCAF.like, 2)
# get cluster means
aggregate(d$CAFs.myCAF.like,by=list(prop_cluster$cluster),FUN=mean)
d$CAFs.myCAF.like_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ CAFs.myCAF.like_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$CAFs.myCAF.like_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


##################### Survival analysis on proportion of Cancer.Basal.SC ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Cancer.Basal.SC), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Cancer.Basal.SC), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Cancer.Basal.SC, 2)
# get cluster means
aggregate(d$Cancer.Basal.SC,by=list(prop_cluster$cluster),FUN=mean)
d$Cancer.Basal.SC_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Cancer.Basal.SC_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Cancer.Basal.SC_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Cancer.Cycling ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Cancer.Cycling), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Cancer.Cycling), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Cancer.Cycling, 2)
# get cluster means
aggregate(d$Cancer.Cycling,by=list(prop_cluster$cluster),FUN=mean)
d$Cancer.Cycling_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Cancer.Cycling_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Cancer.Cycling_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Cancer.Her2.SC ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Cancer.Her2.SC), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Cancer.Her2.SC), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Cancer.Her2.SC, 2)
# get cluster means
aggregate(d$Cancer.Her2.SC,by=list(prop_cluster$cluster),FUN=mean)
d$Cancer.Her2.SC_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Cancer.Her2.SC_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Cancer.Her2.SC_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


##################### Survival analysis on proportion of Cancer.LumA.SC ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Cancer.LumA.SC), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Cancer.LumA.SC), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Cancer.LumA.SC, 2)
# get cluster means
aggregate(d$Cancer.LumA.SC,by=list(prop_cluster$cluster),FUN=mean)
d$Cancer.LumA.SC_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Cancer.LumA.SC_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Cancer.LumA.SC_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Cancer.LumB.SC ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Cancer.LumB.SC), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Cancer.LumB.SC), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Cancer.LumB.SC, 2)
# get cluster means
aggregate(d$Cancer.LumB.SC,by=list(prop_cluster$cluster),FUN=mean)
d$Cancer.LumB.SC_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Cancer.LumB.SC_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Cancer.LumB.SC_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Cycling_Myeloid ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Cycling_Myeloid), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Cycling_Myeloid), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Cycling_Myeloid, 2)
# get cluster means
aggregate(d$Cycling_Myeloid,by=list(prop_cluster$cluster),FUN=mean)
d$Cycling_Myeloid_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Cycling_Myeloid_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Cycling_Myeloid_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Cycling.PVL ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Cycling.PVL), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Cycling.PVL), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Cycling.PVL, 2)
# get cluster means
aggregate(d$Cycling.PVL,by=list(prop_cluster$cluster),FUN=mean)
d$Cycling.PVL_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Cycling.PVL_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Cycling.PVL_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


##################### Survival analysis on proportion of Cycling.T.cells ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Cycling.T.cells), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Cycling.T.cells), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Cycling.T.cells, 2)
# get cluster means
aggregate(d$Cycling.T.cells,by=list(prop_cluster$cluster),FUN=mean)
d$Cycling.T.cells_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Cycling.T.cells_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Cycling.T.cells_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of DCs ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$DCs), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$DCs), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$DCs, 2)
# get cluster means
aggregate(d$DCs,by=list(prop_cluster$cluster),FUN=mean)
d$DCs_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ DCs_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$DCs_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Endothelial.ACKR1 ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Endothelial.ACKR1), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Endothelial.ACKR1), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Endothelial.ACKR1, 2)
# get cluster means
aggregate(d$Endothelial.ACKR1,by=list(prop_cluster$cluster),FUN=mean)
d$Endothelial.ACKR1_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Endothelial.ACKR1_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Endothelial.ACKR1_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


##################### Survival analysis on proportion of Endothelial.CXCL12 ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Endothelial.CXCL12), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Endothelial.CXCL12), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Endothelial.CXCL12, 2)
# get cluster means
aggregate(d$Endothelial.CXCL12,by=list(prop_cluster$cluster),FUN=mean)
d$Endothelial.CXCL12_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Endothelial.CXCL12_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Endothelial.CXCL12_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


##################### Survival analysis on proportion of Endothelial.Lymphatic.LYVE1 ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Endothelial.Lymphatic.LYVE1), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Endothelial.Lymphatic.LYVE1), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Endothelial.Lymphatic.LYVE1, 2)
# get cluster means
aggregate(d$Endothelial.Lymphatic.LYVE1,by=list(prop_cluster$cluster),FUN=mean)
d$Endothelial.Lymphatic.LYVE1_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Endothelial.Lymphatic.LYVE1_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Endothelial.Lymphatic.LYVE1_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


##################### Survival analysis on proportion of Endothelial.RGS5 ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Endothelial.RGS5), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Endothelial.RGS5), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Endothelial.RGS5, 2)
# get cluster means
aggregate(d$Endothelial.RGS5,by=list(prop_cluster$cluster),FUN=mean)
d$Endothelial.RGS5_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Endothelial.RGS5_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Endothelial.RGS5_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Luminal.Progenitors ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Luminal.Progenitors), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Luminal.Progenitors), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Luminal.Progenitors, 2)
# get cluster means
aggregate(d$Luminal.Progenitors,by=list(prop_cluster$cluster),FUN=mean)
d$Luminal.Progenitors_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Luminal.Progenitors_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Luminal.Progenitors_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


##################### Survival analysis on proportion of Macrophage ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Macrophage), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Macrophage), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Macrophage, 2)
# get cluster means
aggregate(d$Macrophage,by=list(prop_cluster$cluster),FUN=mean)
d$Macrophage_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Macrophage_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Macrophage_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Mature.Luminal ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Mature.Luminal), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Mature.Luminal), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Mature.Luminal, 2)
# get cluster means
aggregate(d$Mature.Luminal,by=list(prop_cluster$cluster),FUN=mean)
d$Mature.Luminal_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Mature.Luminal_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Mature.Luminal_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Monocyte ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Monocyte), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Monocyte), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Monocyte, 2)
# get cluster means
aggregate(d$Monocyte,by=list(prop_cluster$cluster),FUN=mean)
d$Monocyte_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Monocyte_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Monocyte_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Myoepithelial ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Myoepithelial), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Myoepithelial), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Myoepithelial, 2)
# get cluster means
aggregate(d$Myoepithelial,by=list(prop_cluster$cluster),FUN=mean)
d$Myoepithelial_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Myoepithelial_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Myoepithelial_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of NK.cells ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$NK.cells), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$NK.cells), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$NK.cells, 2)
# get cluster means
aggregate(d$NK.cells,by=list(prop_cluster$cluster),FUN=mean)
d$NK.cells_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ NK.cells_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$NK.cells_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of NKT.cells ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$NKT.cells), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$NKT.cells), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$NKT.cells, 2)
# get cluster means
aggregate(d$NKT.cells,by=list(prop_cluster$cluster),FUN=mean)
d$NKT.cells_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ NKT.cells_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$NKT.cells_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of Plasmablasts ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Plasmablasts), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Plasmablasts), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Plasmablasts, 2)
# get cluster means
aggregate(d$Plasmablasts,by=list(prop_cluster$cluster),FUN=mean)
d$Plasmablasts_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Plasmablasts_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Plasmablasts_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of PVL.Differentiated ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$PVL.Differentiated), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$PVL.Differentiated), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$PVL.Differentiated, 2)
# get cluster means
aggregate(d$PVL.Differentiated,by=list(prop_cluster$cluster),FUN=mean)
d$PVL.Differentiated_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ PVL.Differentiated_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$PVL.Differentiated_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


##################### Survival analysis on proportion of PVL.Immature ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$PVL.Immature), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$PVL.Immature), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$PVL.Immature, 2)
# get cluster means
aggregate(d$PVL.Immature,by=list(prop_cluster$cluster),FUN=mean)
d$PVL.Immature_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ PVL.Immature_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$PVL.Immature_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of T.cells.CD4. ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$T.cells.CD8.), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$T.cells.CD8.), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$T.cells.CD8., 2)
# get cluster means
aggregate(d$T.cells.CD8.,by=list(prop_cluster$cluster),FUN=mean)
d$T.cells.CD8._kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ T.cells.CD8._kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$T.cells.CD8._kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##################### Survival analysis on proportion of T.cells.CD8. ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$T.cells.CD4.), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$T.cells.CD4.), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$T.cells.CD4., 2)
# get cluster means
aggregate(d$T.cells.CD4.,by=list(prop_cluster$cluster),FUN=mean)
d$T.cells.CD4._kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ T.cells.CD4._kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$T.cells.CD4._kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


### Survival analysis from https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/
library(ggfortify)
# Fit Cox Model
cox <- coxph(Surv(time, vital_status) ~ Basal + Basal_Myoepithelial + Luminal_1_1 + Luminal_1_2 + Luminal_2 + Unclassified +
               B_cell + Fibroblast + Monocyte + T_cell + Topic1 + Topic2 + Topic3 + Topic4 + Topic5 + Topic6 + Topic7 +
               Topic8 + Topic9 + Topic10 , data = d)
summary(cox)
cox_fit <- survfit(cox)
autoplot(cox_fit)

aa_fit <-aareg(Surv(time, vital_status) ~ Basal + Basal_Myoepithelial + Luminal_1_1 + Luminal_1_2 + Luminal_2 + Unclassified +
                 B_cell + Fibroblast + Monocyte + T_cell + Topic1 + Topic2 + Topic3 + Topic4 + Topic5 + Topic6 + Topic7 +
                 Topic8 + Topic9 + Topic10 , data = d)

### GSEA analysis for bulkRS topics
library(fgsea)
library(msigdbr)
library(ggplot2)

# Get genesets for pancreatic tissue from CellMarker Database
cellmarker <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/Human_cell_markers_geneset.txt",sep="#")
cellMarker_gs <- split(x = cellmarker$GeneSymbol, f = cellmarker$TissueCellName)

#Hallmark gene sets
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
h_gene_sets = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

c2_gene_sets = msigdbr(species = "Homo sapiens", category = "C2")
c2_gene_sets = split(x = c2_gene_sets$gene_symbol, f = c2_gene_sets$gs_name)

f <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/TCGA/BRCA/unsup/trainData_JCVB0_nmar_K10_iter100_phi_normalized.csv",header=FALSE)
genes <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/TCGA/BRCA/unsup/genes_bulk_brca",header=FALSE)

df <- cbind(genes,f)

rownames(df)<-df$V1

## Topic1
gseaDat <- df
ranks <- gseaDat$V3
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 

gseaDat <- df
ranks <- gseaDat$V3
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 


## Topic3
gseaDat <- df
ranks <- gseaDat$V5
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 

gseaDat <- df
ranks <- gseaDat$V5
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 

## Topic4
gseaDat <- df
ranks <- gseaDat$V6
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 

gseaDat <- df
ranks <- gseaDat$V6
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 


## Topic7
gseaDat <- df
ranks <- gseaDat$V9
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 

gseaDat <- df
ranks <- gseaDat$V9
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 


## Topic8
gseaDat <- df
ranks <- gseaDat$V10
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 

gseaDat <- df
ranks <- gseaDat$V10
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
