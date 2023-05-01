# Fig5: GTM-decon for phenotype-guided inference
setwd("~/Swapna/WORK/McGill/projects/GTM_decon/manuscript/GB_revision/post_processing/")

library(pheatmap)
library(reshape2)
library(ggplot2)
library(ggpubr)

library(forcats)
library(dplyr)
library(Rmisc)


## Snippet for Converting bulkRS to sparse matrix

#f<-read.csv(file=<path to bulkRNAseq file>,sep="\t",row.names=1)
#f2<-f[,-c(ncol(f))]
#avg_val <- apply(f2,1,quantile,c(0.75))
#f[f<avg_val] <-0
#write.table(f,file="<Path to output file>",sep="\t")


## Fig5a: Phi matrix for Basal vs. ER+ phenotypes

f <- read.csv(file="../data/Fig5_BRCA_Basal_ER_trainData_5topics_phi_normalized.csv",header=FALSE)
genes <-read.csv(file="../data/Fig5_BRCA_Basal_ER_genes.txt",header=FALSE)

df <- cbind(genes,f)
#rownames(df)<-df$V1

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
categories <- c("Basal - Topic1","B - Topic2","B - Topic3","B - Topic4","B - Topic5","ER+ - Topic1","E - Topic2","E - Topic3","E - Topic4","E - Topic5")

## DESEQ2 for ER+ vs TNBC
f <- read.csv(file="../data/Fig5_DEgenes_BRCA_ER_vs_TNBC_together.csv",header=TRUE)
#pheatmap(as.matrix(f[9:18]),cluster_rows = F, cluster_cols = T,color=colorRampPalette(c("orange","blue"))(50),show_rownames=F,cellwidth=20,fontsize = 8)
#log2FC
#pheatmap(f[2],cluster_rows = F, cluster_cols = F,color=colorRampPalette(c("red","white","blue"))(50),cellwidth = 20,show_rownames = F)
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
png("../data/Fig5a.png", units="in", width=7, height=7, res=300)
pheatmap(matrix,annotation_row = annot_row, annotation_colors = ann_colors, cluster_cols=F,cluster_rows=F,color=colorRampPalette(c("white","blue"))(50),show_rownames = F,fontsize = 12,labels_col = categories)
dev.off()

## Fig5b - Deconvolution for held-out dataset
f <-read.csv(file="../data/Fig5_BRCA_Clinical_notes.csv")
subsetf <- read.csv(file="../data/Fig5_BRCA_Basal_ER_sparse_samples.lst",header=FALSE)
f <- f[match(subsetf$V1,f$SampleID),]

#Get the correct order for shuffled set - for clinical data
train_order <- read.csv(file="../data/Fig5_BRCA_Basal_ER_sparse_train.csv",header=FALSE)
test_order <- read.csv(file="../data/Fig5_BRCA_Basal_ER_sparse_test.csv",header=FALSE)

o<-as.integer(test_order$V1)
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

m <- read.csv(file="../data/Fig5_BRCA_Basal_ER_testData_5topics_metagene_normalized.csv",header=FALSE)
row.names(m) <- row.names(clinical)
ordered_matrix <- m[order(m[,2]),]
categories <- c("Basal","ER+")
png("../data/Fig5b.png", units="in", width=7, height=7, res=300)
pheatmap(ordered_matrix[,c(1,2)],annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 12,fontsize_col = 12,cluster_cols = F, cluster_rows = F, labels_col = categories,cellwidth = 40)
dev.off()


## BRCA - histological type - GTM

## Fig5e
f <-read.csv(file="../data/Fig5_BRCA_Clinical_notes.csv")
subsetf <- read.csv(file="../data/Fig5_BRCA_histological_type_sparse_samples.lst",header=FALSE)
f <- f[match(subsetf$V1,f$SampleID),]

#Get the correct order for shuffled set - for clinical data
train_order <- read.csv(file="../data/Fig5_BRCA_histological_type_sparse_train.csv",header=FALSE)
test_order <- read.csv(file="../data/Fig5_BRCA_histological_type_sparse_test.csv",header=FALSE)

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

m <- read.csv(file="../data/Fig5_BRCA_histological_type_testData_5topics_metagene_normalized.csv",header=FALSE)
row.names(m) <- row.names(clinical)
ordered_matrix <- m[order(m[,2]),]
categories <- c("ductal carcinoma","lobular carcinoma")
png("../data/Fig5e.png", units="in", width=7, height=7, res=300)
pheatmap(ordered_matrix[,c(1,2)],annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 10,fontsize_col = 12,cluster_cols = F, cluster_rows = F, labels_col = categories,cellwidth = 40)
dev.off()

## Fig5d: Phi matrix - histological types

f <- read.csv(file="../data/Fig5_BRCA_histological_type_trainData_5topics_phi_normalized.csv",header=FALSE)
genes <-read.csv(file="../data/Fig5_BRCA_histological_type_genes.txt",header=FALSE)

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
f <- read.csv(file="../data/Fig5_DESEQ2_BRCA_type_lobular_vs_ductal_lfc0.58_padj0.1.csv",header=TRUE)
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
png("../data/Fig5d.png", units="in", width=7, height=7, res=300)
pheatmap(matrix,annotation_row = annot_row, annotation_colors = ann_colors, cluster_cols=F,cluster_rows=F,color=colorRampPalette(c("white","blue"))(50),show_rownames = F,fontsize = 12,labels_col = categories)
dev.off()


# GSEA analysis

library(fgsea)
library(msigdbr)

#Hallmark gene sets
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
h_gene_sets = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

f <- read.csv(file="../data/Fig5_BRCA_Basal_ER_trainData_5topics_phi_normalized.csv",header=FALSE)
genes <-read.csv(file="../data/Fig5_Basal_ER_genes.txt",header=FALSE)

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