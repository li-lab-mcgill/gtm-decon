#bulkRS_GTM_BRCA

# GTM for bulkRS based on clinical data

library(pheatmap)
library(reshape2)
library(ggplot2)
library(ggpubr)

library(forcats)
library(dplyr)
library(Rmisc)


## Convert bulkRS to sparse matrix
coding<-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/HGNC_protein_coding_genes.csv",header=FALSE)
coding_genes<-coding$V1
#coding_genes<-append(coding_genes,c("Unnamed:0","histological_type"))
coding_genes<-append(coding_genes,c("Unnamed:0","pathologic_stage"))


f<-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA.pp.TNBC_ER.tab",sep="\t",row.names=1)

#f <- f[intersect(names(f),coding_genes)]
f2<-f[,-c(ncol(f))]
#f2<-f[,-c(20506,20505,20504,20503)] ## For BRCA.pp.TNBC_ER.tab - i dont know why
avg_val <- apply(f2,1,quantile,c(0.75))
#avg_val <- mean(apply(f2,1,quantile,c(0.75)))
#avg_val <- mean(apply(f2,1,quantile,c(0.50)))
f[f<avg_val] <-0
write.table(f,file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA.pp.TNBC_ER.sparse.75p.tab",sep="\t")


## TNBC - ER - GTM
f <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_Clinical_notes_extra.csv")
subsetf <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_TNBC_ER_sparse_75p_samples.lst",header=FALSE)
f <- f[match(subsetf$V1,f$SampleID),]

#Get the correct order for shuffled set - for clinical data
train_order <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_TNBC_ER_sparse_75p_train.csv",header=FALSE)
test_order <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_TNBC_ER_sparse_75p_test.csv",header=FALSE)

o<-as.integer(test_order$V1)
#o<-as.integer(train_order$V1)
o<-o+1
clinical <- f[o,]

clinical$TNBC_ER <- gsub("TNBC", "Basal", clinical$TNBC_ER)
annot_row <- data.frame("Basal_ER"=clinical$TNBC_ER)
row.names(annot_row) <- row.names(clinical)


annot_colors = list(
  #subtype=c(infiltrating_lobular_carcinoma = "#66C2A5", infiltrating_ductal_carcinoma = "#D3D3D3", infiltrating_carcinoma_nos = "#8DA0CB", mucinous_carcinoma = "#E78AC3", medullary_carcinoma = "#A6D854", metaplastic_carcinoma = "#FFD92F", mixed_histology = "#E5C494", other = "#B3B3B3", na = "#FFFFFF"),
  subtype=c(infiltrating_lobular_carcinoma = "#008000", infiltrating_ductal_carcinoma = "#D3D3D3"),
  gender=c(female = "#66A61E", male = "#E6AB02"),
  stage=c(stage_i = "#66C2A5",stage_ii = "#FC8D62",stage_iii = "#8DA0CB",stage_iv = "#E78AC3",stage_x = "#A6D854", na = "#FFFFFF"),
  T_stage=c(t1 = "#FFFFCC",t2 = "#FED976",t3 = "#FD8D3C",t4 = "#E31A1C",tx = "#FFFFFF", na="#FFFFFF"),
  M_stage=c(m0 = "#FFFFFF", m1 ="#969696",mx = "#FFFFFF"),
  N_stage=c(n0 = "#FFFFD9",n1 = "#C7E9B4", n2 = "#41B6C4", n3 ="#225EA8", nx = "#FFFFFF", na="#FFFFFF"),
  days_to_death=c("darkgreen","white"),
  PAM50 = c(Basal = "red", Her2 = "green", LumA = "lightblue", LumB = "darkblue", Normal = "yellow", na = "white"),
  ER_IHC = c(Positive = "darkgreen", Negative = "red", Indeterminate = "grey", Not_Evaluated = "grey", na = "white"),
  Basal_ER = c(ER = "lightpink", Basal = "purple")
  #cell_type=c(Breast = "#0000FF", Immune = "#00FF00", Topic = "#A9A9A9")
)

m <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/TCGA/BRCA/gtm_TNBC_ER_sparse_75p/5topics/testData_trainData_JCVB0_nmar_K10_iter4_metagene_normalized.csv",header=FALSE)
#m <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_GTM/TCGA/BRCA/gtm_TNBC_ER_sparse_75p/5topics/pp/testData_trainData_JCVB0_nmar_K10_iter4_metaphe.csv",header=FALSE)
row.names(m) <- row.names(clinical)
ordered_matrix <- m[order(m[,2]),]
categories <- c("Basal","ER+")
pheatmap(ordered_matrix[,c(1,2)],annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 12,fontsize_col = 12,cluster_cols = F, cluster_rows = F, labels_col = categories)
categories <- c("Basal - T1", "T2", "T3", "T4", "T5","ER+ - T1", "T2", "T3", "T4", "T5")
pheatmap(m[1:10],annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 12,fontsize_col = 12,cluster_cols = F, labels_col = categories)


## TNBC - ER - DESEQ2
f <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_Clinical_notes_extra.csv")
subsetf <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_TNBC_ER_DESEQ_samples.lst",header=FALSE)
f <- f[match(subsetf$V1,f$SampleID),]

#Get the correct order for shuffled set - for clinical data
train_order <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_TNBC_ER_DESEQ_train.csv",header=FALSE)
test_order <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_TNBC_ER_DESEQ_test.csv",header=FALSE)

o<-as.integer(test_order$V1)
#o<-as.integer(train_order$V1)
o<-o+1
clinical <- f[o,]

annot_row <- data.frame("TNBC_ER"=clinical$TNBC_ER)
row.names(annot_row) <- row.names(clinical)


annot_colors = list(
  #subtype=c(infiltrating_lobular_carcinoma = "#66C2A5", infiltrating_ductal_carcinoma = "#D3D3D3", infiltrating_carcinoma_nos = "#8DA0CB", mucinous_carcinoma = "#E78AC3", medullary_carcinoma = "#A6D854", metaplastic_carcinoma = "#FFD92F", mixed_histology = "#E5C494", other = "#B3B3B3", na = "#FFFFFF"),
  subtype=c(infiltrating_lobular_carcinoma = "#008000", infiltrating_ductal_carcinoma = "#D3D3D3"),
  gender=c(female = "#66A61E", male = "#E6AB02"),
  stage=c(stage_i = "#66C2A5",stage_ii = "#FC8D62",stage_iii = "#8DA0CB",stage_iv = "#E78AC3",stage_x = "#A6D854", na = "#FFFFFF"),
  T_stage=c(t1 = "#FFFFCC",t2 = "#FED976",t3 = "#FD8D3C",t4 = "#E31A1C",tx = "#FFFFFF", na="#FFFFFF"),
  M_stage=c(m0 = "#FFFFFF", m1 ="#969696",mx = "#FFFFFF"),
  N_stage=c(n0 = "#FFFFD9",n1 = "#C7E9B4", n2 = "#41B6C4", n3 ="#225EA8", nx = "#FFFFFF", na="#FFFFFF"),
  days_to_death=c("darkgreen","white"),
  PAM50 = c(Basal = "red", Her2 = "green", LumA = "lightblue", LumB = "darkblue", Normal = "yellow", na = "white"),
  ER_IHC = c(Positive = "darkgreen", Negative = "red", Indeterminate = "grey", Not_Evaluated = "grey", na = "white"),
  TNBC_ER = c(ER = "lightpink", TNBC = "purple")
  #cell_type=c(Breast = "#0000FF", Immune = "#00FF00", Topic = "#A9A9A9")
)

m <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_GTM/TCGA/BRCA/gtm_TNBC_ER_DESEQ2//5topics/pp/testData_trainData_JCVB0_nmar_K10_iter4_metaphe_normalized.csv",header=FALSE)
m <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_GTM/TCGA/BRCA/gtm_TNBC_ER_DESEQ2/5topics/pp/testData_trainData_JCVB0_nmar_K10_iter4_metaphe.csv",header=FALSE)
row.names(m) <- row.names(clinical)
ordered_matrix <- m[order(m[,2]),]
categories <- c("Basal","ER+")
pheatmap(ordered_matrix[,c(1,2)],annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 12,fontsize_col = 12,cluster_cols = F, cluster_rows = F, labels_col = categories)
categories <- c("Basal - T1", "B - T2", "B - T3", "B - T4", "B - T5","ER+ - T1", "E - T2", "E - T3", "E - T4", "E - T5")
pheatmap(m[1:10],annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 12,fontsize_col = 12,cluster_cols = F, labels_col = categories)



## BRCA - histological type - GTM
f <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_Clinical_notes_extra.csv")
subsetf <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_type_2sets_sparse_75p_samples.lst",header=FALSE)
f <- f[match(subsetf$V1,f$SampleID),]

#Get the correct order for shuffled set - for clinical data
train_order <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_type_2sets_sparse_75p_train.csv",header=FALSE)
test_order <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_type_2sets_sparse_75p_test.csv",header=FALSE)

o<-as.integer(test_order$V1)
#o<-as.integer(train_order$V1)
o<-o+1
clinical <- f[o,]

annot_row <- data.frame("subtype"=clinical$histological_type)
row.names(annot_row) <- row.names(clinical)


annot_colors = list(
  #subtype=c(infiltrating_lobular_carcinoma = "#66C2A5", infiltrating_ductal_carcinoma = "#D3D3D3", infiltrating_carcinoma_nos = "#8DA0CB", mucinous_carcinoma = "#E78AC3", medullary_carcinoma = "#A6D854", metaplastic_carcinoma = "#FFD92F", mixed_histology = "#E5C494", other = "#B3B3B3", na = "#FFFFFF"),
  subtype=c(infiltrating_lobular_carcinoma = "#008000", infiltrating_ductal_carcinoma = "#D3D3D3"),
  gender=c(female = "#66A61E", male = "#E6AB02"),
  stage=c(stage_i = "#66C2A5",stage_ii = "#FC8D62",stage_iii = "#8DA0CB",stage_iv = "#E78AC3",stage_x = "#A6D854", na = "#FFFFFF"),
  T_stage=c(t1 = "#FFFFCC",t2 = "#FED976",t3 = "#FD8D3C",t4 = "#E31A1C",tx = "#FFFFFF", na="#FFFFFF"),
  M_stage=c(m0 = "#FFFFFF", m1 ="#969696",mx = "#FFFFFF"),
  N_stage=c(n0 = "#FFFFD9",n1 = "#C7E9B4", n2 = "#41B6C4", n3 ="#225EA8", nx = "#FFFFFF", na="#FFFFFF"),
  days_to_death=c("darkgreen","white"),
  PAM50 = c(Basal = "red", Her2 = "green", LumA = "lightblue", LumB = "darkblue", Normal = "yellow", na = "white"),
  ER_IHC = c(Positive = "darkgreen", Negative = "red", Indeterminate = "grey", Not_Evaluated = "grey", na = "white"),
  Basal_ER = c(ER = "grey", Basal = "purple", na = "white")
  #cell_type=c(Breast = "#0000FF", Immune = "#00FF00", Topic = "#A9A9A9")
)

row.names(annot_row)<-row.names(clinical)

m <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_GTM/TCGA/BRCA/gtm_type_2sets_sparse//5topics/pp/trainData_trainData_JCVB0_nmar_K10_iter19_metaphe.csv",header=FALSE)
m <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_GTM/TCGA/BRCA/gtm_type_2sets_sparse//5topics/pp/trainData_trainData_JCVB0_nmar_K10_iter19_metaphe_normalized.csv",header=FALSE)
m <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/TCGA/BRCA/gtm_histological_type_sparse_75p//5topics/testData_trainData_JCVB0_nmar_K10_iter19_metagene_normalized.csv",header=FALSE)
row.names(m) <- row.names(clinical)
ordered_matrix <- m[order(m[,2]),]
categories <- c("ductal carcinoma","lobular carcinoma")
pheatmap(ordered_matrix[,c(1,2)],annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 10,fontsize_col = 12,cluster_cols = F, cluster_rows = F, labels_col = categories)
categories <- c("ductal carcinoma - T1", "T2", "T3", "T4", "T5","lobular carcinoma - T1", "T2", "T3", "T4", "T5")
pheatmap(m[1:10],annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 8,fontsize_col = 10,cluster_cols = F, labels_col = categories)


## BRCA - histological type - DESEQ2
f <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_Clinical_notes_extra.csv")
subsetf <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_type_DESEQ_samples.lst",header=FALSE)
f <- f[match(subsetf$V1,f$SampleID),]

#Get the correct order for shuffled set - for clinical data
train_order <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_type_DESEQ_train.csv",header=FALSE)
test_order <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/BRCA_type_DESEQ_test.csv",header=FALSE)

#o<-as.integer(test_order$V1)
o<-as.integer(train_order$V1)
o<-o+1
clinical <- f[o,]

annot_row <- data.frame("subtype"=clinical$histological_type)
row.names(annot_row) <- row.names(clinical)


annot_colors = list(
  #subtype=c(infiltrating_lobular_carcinoma = "#66C2A5", infiltrating_ductal_carcinoma = "#D3D3D3", infiltrating_carcinoma_nos = "#8DA0CB", mucinous_carcinoma = "#E78AC3", medullary_carcinoma = "#A6D854", metaplastic_carcinoma = "#FFD92F", mixed_histology = "#E5C494", other = "#B3B3B3", na = "#FFFFFF"),
  subtype=c(infiltrating_lobular_carcinoma = "#008000", infiltrating_ductal_carcinoma = "#D3D3D3"),
  gender=c(female = "#66A61E", male = "#E6AB02"),
  stage=c(stage_i = "#66C2A5",stage_ii = "#FC8D62",stage_iii = "#8DA0CB",stage_iv = "#E78AC3",stage_x = "#A6D854", na = "#FFFFFF"),
  T_stage=c(t1 = "#FFFFCC",t2 = "#FED976",t3 = "#FD8D3C",t4 = "#E31A1C",tx = "#FFFFFF", na="#FFFFFF"),
  M_stage=c(m0 = "#FFFFFF", m1 ="#969696",mx = "#FFFFFF"),
  N_stage=c(n0 = "#FFFFD9",n1 = "#C7E9B4", n2 = "#41B6C4", n3 ="#225EA8", nx = "#FFFFFF", na="#FFFFFF"),
  days_to_death=c("darkgreen","white"),
  PAM50 = c(Basal = "red", Her2 = "green", LumA = "lightblue", LumB = "darkblue", Normal = "yellow", na = "white"),
  ER_IHC = c(Positive = "darkgreen", Negative = "red", Indeterminate = "grey", Not_Evaluated = "grey", na = "white"),
  TNBC_ER = c(ER = "grey", TNBC = "purple", na = "white")
  #cell_type=c(Breast = "#0000FF", Immune = "#00FF00", Topic = "#A9A9A9")
)

row.names(annot_row)<-row.names(clinical)

m <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_GTM/TCGA/BRCA/gtm_histological_type_DESEQ2///5topics/pp/trainData_trainData_JCVB0_nmar_K10_iter19_metaphe.csv",header=FALSE)
#m <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_GTM/TCGA/BRCA/gtm_histological_type_DESEQ2//5topics/pp/testData_trainData_JCVB0_nmar_K10_iter19_metaphe_normalized.csv",header=FALSE)
row.names(m) <- row.names(clinical)
ordered_matrix <- m[order(m[,2]),]
categories <- c("ductal carcinoma","lobular carcinoma")
pheatmap(ordered_matrix[,c(1,2)],annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 10,fontsize_col = 12,cluster_cols = F, cluster_rows = F, labels_col = categories)
categories <- c("ductal carcinoma - T1", "T2", "T3", "T4", "T5","lobular carcinoma - T1", "T2", "T3", "T4", "T5")
pheatmap(m[1:10],annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 8,fontsize_col = 10,cluster_cols = F, labels_col = categories)



## Phi matrix - TNBC vs. ER

f <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/TCGA/BRCA/gtm_TNBC_ER_sparse_75p/5topics/trainData_JCVB0_nmar_K10_iter4_phi_normalized.csv",header=FALSE)
genes <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/TCGA/BRCA/gtm_TNBC_ER_sparse_75p/5topics/genes.txt",header=FALSE)

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
#subtypes <- c("infiltrating_ductal_carcinoma", "infiltrating_lobular_carcinoma", "medullary_carcinoma", "metaplastic_carcinoma", "mixed_histology", "mucinous_carcinoma", "other")
#subtypes <- c("ductal carcinoma - T1","ductal carcinoma - T2","ductal carcinoma - T3","ductal carcinoma - T4","lobular carcinoma - T1", "lobular carcinoma - T2", "lobular carcinoma - T3", "lobular carcinoma - T4")
#subtypes <- c("ductal carcinoma - T1","T2","T3","T4","T5","lobular carcinoma - T1", "T2", "T3", "T4", "T5")
#subtypes <- c("ductal carcinoma - T1","ductal carcinoma - T2","ductal carcinoma - T3","ductal carcinoma - T4","ductal carcinoma - T5","lobular carcinoma - T1", "lobular carcinoma - T2", "lobular carcinoma - T3", "lobular carcinoma - T4", "lobular carcinoma - T5")
categories <- c("Basal - Topic1","B - Topic2","B - Topic3","B - Topic4","B - Topic5","ER+ - Topic1","E - Topic2","E - Topic3","E - Topic4","E - Topic5")
#categories <- c("ductal carcinoma - T1","T2","T3","T4","T5","lobular carcinoma - T1","T2","T3","T4","T5")
## Marker genes from PangloDB
#mammary_epithelial_markers <- c("WNT5B","PRLR","CLDN4","CSN3","CSN1S1","RARRES1","KRT17","NNMT","ALDH1A3","SGMS2","PIEZO1","ANO1","KRT5","NEDD9","FOXP4","LSR","LIF","KRT19","KRT7","STAT3","KRT8","KRT18","BTN1A1","BTN2A2","BTN3A1","MMP14","TOX3","CDH1","FFAR1","PYGO2","SERPINB5","RUNX2","CSN2","CITED1","IRF6","SULT1E1","KRT14","KDM3A","SLC30A2","PADI2","CLDN1","ABCG2","PIP")
#luminal_epithelial_markers <- c("AR","PIP","SERPINA1","ATP7B","HIF1A","FGFR4","DDR1","CEBPD","PGR","KRT8","UXT","PTH1R","CD9","AQP3","ATP2C2","WNT5A","SLC12A2","ESR1","AQP5","RUNX1","CDH1","MUC1","ANPEP","KLK3","KRT23","FGG","KRT18","ANKRD30A","SLPI","PROM1","KRT19","SYTL2","CD74","AGR2","LTF","SAA2","FGFR2","SERPINB4","SERPINB3","WFDC2","LCN2","BTG1","CLDN4","ANXA1","HMGA1","STC2","AREG","TNFSF10")

#both_markers <- mammary_epithelial_markers
#both_markers <- append(both_markers, luminal_epithelial_markers)
#a <-topn_df$V1
#annot_list_mammary <- a %in% mammary_epithelial_markers
#a <-topn_df$V1
#annot_list_luminal <- a %in% luminal_epithelial_markers

#annot_row = data.frame(Mammary_markers = factor(annot_list_mammary),
#                       Luminal_markers = factor(annot_list_luminal))
#ann_colors = list(Mammary_markers=c("TRUE" = "red", "FALSE" = "white", "NA" = "lightgrey"),
#                  Luminal_markers = c("TRUE" = "orange", "FALSE" = "white", "NA" = "lightgrey"))

## Most different genes
df_diff <- ((f$V3+f$V4+f$V5+f$V6+f$V7)/5) - ((f$V8+f$V9+f$V10+f$V11+f$V12)/5)
abs_rowsums <- abs(df_diff)
ndf<-cbind(genes,abs_rowsums)
top_diff_genes <- ndf %>% arrange(desc(ndf[,2])) %>% slice(1:1000)

## DESEQ2 for ER+ vs TNBC
f <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/DEgenes_BRCA_ER_vs_TNBC_together.csv",header=TRUE)
#pheatmap(as.matrix(f[9:18]),cluster_rows = F, cluster_cols = T,color=colorRampPalette(c("orange","blue"))(50),show_rownames=F,cellwidth=20,fontsize = 8)
#log2FC
pheatmap(f[2],cluster_rows = F, cluster_cols = F,color=colorRampPalette(c("red","white","blue"))(50),cellwidth = 20,show_rownames = F)
deseq_up_genes <- f[f$log2FoldChange>0,]$Genes
deseq_down_genes <- f[f$log2FoldChange<0,]$Genes


up_genes <- topn_df$V1 %in% deseq_up_genes
down_genes <- topn_df$V1 %in% deseq_down_genes

annot_row = data.frame(DESEQ_ER_Up = factor(up_genes),
                       DESEQ_Basal_Up = factor(down_genes))
ann_colors = list(DESEQ_ER_Up=c("TRUE" = "lightpink", "FALSE" = "white"),
                  DESEQ_Basal_Up = c("TRUE" = "purple", "FALSE" = "white"))


matrix <- as.data.frame(select(topn_df,2:ncol(topn_df)))
rownames(matrix) <- rownames(topn_df)
rownames(annot_row) <- rownames(topn_df)

max_val <- max(matrix)
matrix <- matrix/max_val
#pheatmap(matrix,cluster_cols=F,cluster_rows=F,color=colorRampPalette(c("white","blue"))(50),show_rownames = F,fontsize = 12,labels_col = categories)
pheatmap(matrix,annotation_row = annot_row, annotation_colors = ann_colors, cluster_cols=F,cluster_rows=F,color=colorRampPalette(c("white","blue"))(50),show_rownames = F,fontsize = 12,labels_col = categories)

##Common with DESEQ and top1000 differences
common <- f$Genes %in% top_diff_genes$V1
length(common[common == TRUE]) ## 779

library(fgsea)
library(msigdbr)

#Hallmark gene sets
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")

##Enrichment analysis for top1000 DE genes
## Over-representation analysis
library(clusterProfiler)

genes_of_interest <- top_diff_genes$V1
background_genes <- df$V1
ORA_DE <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = select(h_gene_sets,gs_name,gene_symbol))

## Phi matrix - histological types

f <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/TCGA/BRCA/gtm_histological_type_sparse_75p/5topics/trainData_JCVB0_nmar_K10_iter19_phi_normalized.csv",header=FALSE)
genes <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/TCGA/BRCA/gtm_histological_type_sparse_75p/5topics/genes.txt",header=FALSE)

df <- cbind(genes,f)
#rownames(df)<-df$V1

#Drop the first second and third column - they dont correspond to phi values
df<-df[-(2:3)]

#Get top50 genes per column
datalist = list()
topn_df <- data.frame()
for(i in 2:ncol(df)) {       # for-loop over columns
  subset_df <- df %>% arrange(desc(df[,i])) %>% slice(1:20)
  datalist[[i]]<-subset_df
}
topn_df <- do.call(rbind,datalist)

categories <- c("ductal carcinoma - Topic1","D - Topic2","D - Topic3","D - Topic4","D - Topic5","lobular carcinoma - Topic1","L - Topic2","L - Topic3","L - Topic4","L - Topic5")

## DESEQ2 for histological subtypes - lobular vs. ductal carcinoma
f <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts/figures/DESEQ2_BRCA_type_lobular_vs_ductal_lfc0.58_padj0.1.csv",header=TRUE)
pheatmap(f[3],cluster_rows = F, cluster_cols = F,color=colorRampPalette(c("red","white","blue"))(50),cellwidth = 20,show_rownames = F)
deseq_up_genes <- f[f$log2FoldChange>0,]$Genes
deseq_down_genes <- f[f$log2FoldChange<0,]$Genes

up_genes <- topn_df$V1 %in% deseq_up_genes
down_genes <- topn_df$V1 %in% deseq_down_genes

annot_row = data.frame(DESEQ_Lobular_Up = factor(up_genes),
                       DESEQ_Ductal_Up = factor(down_genes))
ann_colors = list(DESEQ_Lobular_Up=c("TRUE" = "darkgreen", "FALSE" = "white"),
                  DESEQ_Ductal_Up = c("TRUE" = "grey", "FALSE" = "white"))


matrix <- as.data.frame(select(topn_df,2:ncol(topn_df)))
rownames(matrix) <- rownames(topn_df)
rownames(annot_row) <- rownames(topn_df)

max_val <- max(matrix)
matrix <- matrix/max_val
#pheatmap(matrix,cluster_cols=F,cluster_rows=F,color=colorRampPalette(c("white","blue"))(50),show_rownames = F,fontsize = 12,labels_col = categories)
pheatmap(matrix,annotation_row = annot_row, annotation_colors = ann_colors, cluster_cols=F,cluster_rows=F,color=colorRampPalette(c("white","blue"))(50),show_rownames = F,fontsize = 12,labels_col = categories)



# GSEA analysis

library(fgsea)
library(msigdbr)

#Hallmark gene sets
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
h_gene_sets = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

c2_gene_sets = msigdbr(species = "Homo sapiens", category = "C2")
c2_gene_sets = split(x = c2_gene_sets$gene_symbol, f = c2_gene_sets$gs_name)

c5_gene_sets = msigdbr(species = "Homo sapiens", category = "C5")
c5_gene_sets = split(x = c5_gene_sets$gene_symbol, f = c5_gene_sets$gs_name)

c6_gene_sets = msigdbr(species = "Homo sapiens", category = "C6")
c6_gene_sets = split(x = c6_gene_sets$gene_symbol, f = c6_gene_sets$gs_name)

f <- read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/TCGA/BRCA/gtm_TNBC_ER_sparse_75p/5topics/trainData_JCVB0_nmar_K10_iter4_phi_normalized.csv",header=FALSE)
genes <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/TCGA/BRCA/gtm_TNBC_ER_sparse_75p/5topics/genes.txt",header=FALSE)

df <- cbind(genes,f)
rownames(df)<-df$V1

gseaDat <- df
ranks <- gseaDat$V3
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]], ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 

ranks <- gseaDat$V4
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]], ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
plotEnrichment(h_gene_sets[["HALLMARK_MTORC1_SIGNALING"]], ranks) + labs(title="HALLMARK_MTORC1_SIGNALING") + theme(text = element_text(size = 20)) 

ranks <- gseaDat$V5
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]], ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 

ranks <- gseaDat$V6
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]], ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
plotEnrichment(h_gene_sets[["HALLMARK_MTORC1_SIGNALING"]], ranks) + labs(title="HALLMARK_MTORC1_SIGNALING") + theme(text = element_text(size = 20)) 

ranks <- gseaDat$V7
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]], ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
plotEnrichment(h_gene_sets[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]], ranks) + labs(title="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") + theme(text = element_text(size = 20)) 

ranks <- gseaDat$V8
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]], ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
plotEnrichment(h_gene_sets[["HALLMARK_P53_PATHWAY"]], ranks) + labs(title="HALLMARK_P53_PATHWAY") + theme(text = element_text(size = 20)) 

ranks <- gseaDat$V9
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]], ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
plotEnrichment(h_gene_sets[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], ranks) + labs(title="HALLMARK_OXIDATIVE_PHOSPHORYLATION") + theme(text = element_text(size = 20)) 

ranks <- gseaDat$V10
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]], ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
plotEnrichment(h_gene_sets[["HALLMARK_MITOTIC_SPINDLE"]], ranks) + labs(title="HALLMARK_MITOTIC_SPINDLE") + theme(text = element_text(size = 20)) 

ranks <- gseaDat$V11
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]], ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 18)) 

ranks <- gseaDat$V12
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]], ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
