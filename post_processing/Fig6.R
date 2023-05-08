setwd("~/Swapna/WORK/McGill/projects/GTM_decon/manuscript/GB_revision/post_processing/")

## Cell-type specific DEgenes

library(pheatmap)
library(dplyr)
library(gtools)

### Phi matrix
f <- read.csv(file="../data/Fig6_DEgenes_trainData_phi_normalized.csv",header=FALSE)
genes <-read.csv(file="../data/Fig6_genes.txt",header=FALSE)

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

matrix <- as.data.frame(select(topn_df,2:ncol(topn_df)))
rownames(matrix) <- rownames(topn_df)
matrix <- matrix/max(matrix)
colnames(matrix)<-c("Basal","Basal_myoepithelial", "Luminal_1_1", "Luminal_1_2", "Luminal_2","Basal","Basal_myoepithelial", "Luminal_1_1", "Luminal_1_2", "Luminal_2")
png("../data/Fig6a.png", units="in", width=7, height=7, res=300)
pheatmap(matrix,cluster_cols=F,cluster_rows=F,color=colorRampPalette(c("white","blue"))(50),show_rownames = F,fontsize = 12)
dev.off()

# GSEA analysis

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


#Holds differences between the phenotypes for each celltype
celltype1 <- df$V3-df$V8
celltype2 <- df$V4-df$V9
celltype3 <- df$V5-df$V10
celltype4 <- df$V6-df$V11
celltype5 <- df$V7-df$V12

celltype6 <- df$V8-df$V3
celltype7 <- df$V9-df$V4
celltype8 <- df$V10-df$V5
celltype9 <- df$V11-df$V6
celltype10 <- df$V12-df$V7


df_diff <- data.frame("Genes" = genes,
                      "B_Basal" = celltype1,
                      "B_Basal_myoepithelial" = celltype2,
                      "B_Luminal_1_1" = celltype3,
                      "B_Luminal_1_2" = celltype4,
                      "B_Luminal_2" = celltype5,
                      "E_Basal" = celltype6,
                      "E_Basal_myoepithelial" = celltype7,
                      "E_Luminal_1_1" = celltype8,
                      "E_Luminal_1_2" = celltype9,
                      "E_Luminal_2" = celltype10)

##Get top20 different genes in cell-type difference 
datalist = list()
topn_df <- data.frame()
for(i in 2:ncol(df_diff)) {       # for-loop over columns
  subset_df <- df_diff %>% arrange(desc(df_diff[,i])) %>% slice(1:20)
  datalist[[i]]<-subset_df
}
topn_df <- do.call(rbind,datalist)

matrix <- as.data.frame(select(topn_df,2:ncol(topn_df)))
rownames(matrix) <- rownames(topn_df)

matrix <- matrix/max(matrix)
png("../data/Fig6b.png", units="in", width=7, height=7, res=300)
pheatmap(matrix,cluster_cols=F,cluster_rows=F,color=colorRampPalette(c("red","white","blue"))(50),show_rownames = F,fontsize = 12)
dev.off()

## Accuracy in deconvolving train and test sets
## TNBC - ER - DESEQ2
f <-read.csv(file="../data/Fig5_BRCA_Clinical_notes.csv")
subsetf <- read.csv(file="../data/Fig5_BRCA_Basal_ER_sparse_samples.lst",header=FALSE)
f <- f[match(subsetf$V1,f$SampleID),]

#Get the correct order for shuffled set - for clinical data
train_order <- read.csv(file="../data/Fig5_BRCA_Basal_ER_sparse_train.csv",header=FALSE)
test_order <- read.csv(file="../data/Fig5_BRCA_Basal_ER_sparse_test.csv",header=FALSE)

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

m <- read.csv(file="../data/Fig6_testData_metagene_normalized.csv",header=FALSE)
row.names(m) <- row.names(clinical)
ordered_matrix <- m[order(m[,2]),]
categories <- c("Basal","ER+")
png("../data/Fig6c_upperpanel.png", units="in", width=7, height=7, res=300)
pheatmap(ordered_matrix[,c(1,2)],annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 20,fontsize_col = 20,cluster_cols = F, cluster_rows = F, labels_col = categories, cellwidth = 40)
dev.off()

m <- read.csv(file="../data/Fig6_testData_metagene.csv",header=FALSE)
row.names(m) <- row.names(clinical)
categories <- c("Basal - Basal", "B - Basal_Myo", "B - Luminal_1_1", "B - Luminal_1_2", "B - Luminal_2","ER+ - Basal", "E - Basal_Myo", "E - Luminal_1_1", "E - Luminal_1_2", "E - Luminal_2")
png("../data/Fig6c_lowerpanel.png", units="in", width=7, height=7, res=300)
pheatmap(m[1:10],annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 12,fontsize_col = 12,cluster_cols = F, labels_col = categories)
dev.off()

# Commented out: Takes long to run - output file provided below to continue with next steps
# ## Permutation test for calculating p-values
# permutation.test <- function(test_vector){
#   pvals=c()
#   for(j in 1:nrow(test_vector)){
#     pvals[j] = 0
#   }
#   
#   diff_tv <- abs(test_vector$Phenotype1 - test_vector$Phenotype2)
#   
#   num_permutations = 100000
#   for(i in 1:num_permutations){
#     df2 <- test_vector
#     n   <- nrow(test_vector)
#     
#     df2[] <- lapply(test_vector,function(x) x[sample(n)] )
#     diff_rv <- abs(df2$Phenotype1 - df2$Phenotype2)
#     
#     for(j in 1:length(diff_tv)){
#       if(diff_rv[j] > diff_tv[j]){
#         pvals[j] = pvals[j] + 1
#       }
#     }
#   }
#   
#   for(j in 1:nrow(test_vector)){
#     pvals[j] = pvals[j]/num_permutations
#   }
#   
#   return(pvals);
# }
# 
# 
# ## For Basal - EvsT 
# test_vector <- data.frame(Phenotype1=df$V8,Phenotype2=df$V3)
# pvals_Basal <- permutation.test(test_vector)
# length(pvals_Basal[pvals_Basal < 0.05])
# 
# ## For Basal_myoepithelial - EvsT 
# test_vector <- data.frame(Phenotype1=df$V9,Phenotype2=df$V4)
# pvals_Basal_myo <- permutation.test(test_vector)
# length(pvals_Basal_myo[pvals_Basal_myo < 0.05])
# 
# ## For Luminal_1_1 - EvsT 
# test_vector <- data.frame(Phenotype1=df$V10,Phenotype2=df$V5)
# pvals_L11 <- permutation.test(test_vector)
# length(pvals_L11[pvals_L11 < 0.05])
# 
# ## For Luminal_1_2 - EvsT 
# test_vector <- data.frame(Phenotype1=df$V11,Phenotype2=df$V6)
# pvals_L12 <- permutation.test(test_vector)
# length(pvals_L12[pvals_L12 < 0.05])
# 
# ## For Luminal_2 - EvsT 
# test_vector <- data.frame(Phenotype1=df$V12,Phenotype2=df$V7)
# pvals_L2 <- permutation.test(test_vector)
# length(pvals_L2[pvals_L2 < 0.05])
# 
# all_pvals <- data.frame(Genes=genes,pvals_Basal=pvals_Basal,pvals_Basal_myo=pvals_Basal_myo,pvals_L11=pvals_L11,pvals_L12=pvals_L12,pvals_L2=pvals_L2)
# all_pvals_lt0.05 <- all_pvals %>% filter_all(any_vars(. < c(.05)))

all_pvals_lt0.05 <- read.csv(file="../data/Fig6_all_pvals_lt0.05.csv")
mat_de <- df[df$V1 %in% all_pvals_lt0.05$V1,]

# ## Over-representation analysis
# library(clusterProfiler)
# 
# genes_pvals_Basal <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_Basal < 0.05,]
# genes_of_interest <- genes_pvals_Basal$V1
# background_genes <- all_pvals$V1
# ORA_Basal <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 
# genes_pvals_Basal_myo <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_Basal_myo < 0.05,]
# genes_of_interest <- genes_pvals_Basal_myo$V1
# background_genes <- all_pvals$V1
# ORA_Basal_myo <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 
# genes_pvals_L11 <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_L11 < 0.05,]
# genes_of_interest <- genes_pvals_L11$V1
# background_genes <- all_pvals$V1
# ORA_L11 <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 
# genes_pvals_L12 <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_L12 < 0.05,]
# genes_of_interest <- genes_pvals_L12$V1
# background_genes <- all_pvals$V1
# ORA_L12 <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 
# genes_pvals_L2 <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_L2 < 0.05,]
# genes_of_interest <- genes_pvals_L2$V1
# background_genes <- all_pvals$V1
# ORA_L2 <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 
# ## COnsidering only upregulated genes with pval < 0.05
# genes_pvals_Basal <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_Basal < 0.05,]
# genes_of_interest <- genes_pvals_Basal[genes_pvals_Basal$V1 %in% (df_diff[df_diff$T_Basal>0,]$V1),]$V1
# background_genes <- all_pvals$V1
# ORA_Basal <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 
# genes_pvals_Basal_myo <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_Basal_myo < 0.05,]
# genes_of_interest <- genes_pvals_Basal_myo[genes_pvals_Basal_myo$V1 %in% (df_diff[df_diff$T_Basal_myoepithelial>0,]$V1),]$V1
# background_genes <- all_pvals$V1
# ORA_Basal_myo <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 
# genes_pvals_L11 <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_L11 < 0.05,]
# genes_of_interest <- genes_pvals_L11[genes_pvals_L11$V1 %in% (df_diff[df_diff$T_Luminal_1_1>0,]$V1),]$V1
# background_genes <- all_pvals$V1
# ORA_L11 <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 
# genes_pvals_L12 <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_L12 < 0.05,]
# genes_of_interest <- genes_pvals_L12[genes_pvals_L12$V1 %in% (df_diff[df_diff$T_Luminal_1_2>0,]$V1),]$V1
# background_genes <- all_pvals$V1
# ORA_L12 <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 
# genes_pvals_L2 <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_L2 < 0.05,]
# genes_of_interest <- genes_pvals_L2[genes_pvals_L2$V1 %in% (df_diff[df_diff$T_Luminal_2>0,]$V1),]$V1
# background_genes <- all_pvals$V1
# ORA_L2 <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 
# ## COnsidering only downregulated genes with pval < 0.05
# genes_pvals_Basal <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_Basal < 0.05,]
# genes_of_interest <- genes_pvals_Basal[genes_pvals_Basal$V1 %in% (df_diff[df_diff$T_Basal<0,]$V1),]$V1
# background_genes <- all_pvals$V1
# ORA_Basal <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 
# genes_pvals_Basal_myo <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_Basal_myo < 0.05,]
# genes_of_interest <- genes_pvals_Basal_myo[genes_pvals_Basal_myo$V1 %in% (df_diff[df_diff$T_Basal_myoepithelial<0,]$V1),]$V1
# background_genes <- all_pvals$V1
# ORA_Basal_myo <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 
# genes_pvals_L11 <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_L11 < 0.05,]
# genes_of_interest <- genes_pvals_L11[genes_pvals_L11$V1 %in% (df_diff[df_diff$T_Luminal_1_1<0,]$V1),]$V1
# background_genes <- all_pvals$V1
# ORA_L11 <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 
# genes_pvals_L12 <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_L12 < 0.05,]
# genes_of_interest <- genes_pvals_L12[genes_pvals_L12$V1 %in% (df_diff[df_diff$T_Luminal_1_2<0,]$V1),]$V1
# background_genes <- all_pvals$V1
# ORA_L12 <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 
# genes_pvals_L2 <- all_pvals_lt0.05[all_pvals_lt0.05$pvals_L2 < 0.05,]
# genes_of_interest <- genes_pvals_L2[genes_pvals_L2$V1 %in% (df_diff[df_diff$T_Luminal_2<0,]$V1),]$V1
# background_genes <- all_pvals$V1
# ORA_L2 <- enricher(genes_of_interest, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_genes, TERM2GENE = dplyr::select(h_gene_sets,gs_name,gene_symbol))
# 


#DEgenes from DESEQ2 for comparison
deseq2 <- read.csv(file="../data/Fig5_DEgenes_BRCA_ER_vs_TNBC_together.csv",header=TRUE)
deseq2_annot <- all_pvals_lt0.05$V1 %in% deseq2$Genes
annot_row = data.frame(DESEQ2 = factor(deseq2_annot))
ann_colors = list(DESEQ2=c("TRUE" = "red", "FALSE" = "white"))
rownames(annot_row)<-rownames(mat_de)
colnames(mat_de)<-c("Genes","Basal","Basal_myoepithelial", "Luminal_1_1", "Luminal_1_2", "Luminal_2","Basal","Basal_myoepithelial", "Luminal_1_1", "Luminal_1_2", "Luminal_2")

mat_de_phi_diff <- df_diff[df_diff$V1 %in% all_pvals_lt0.05$V1,]
colnames(mat_de_phi_diff)<-c("Genes","Basal","Basal_myoepithelial", "Luminal_1_1", "Luminal_1_2", "Luminal_2","Basal","Basal_myoepithelial", "Luminal_1_1", "Luminal_1_2", "Luminal_2")

##Heatmaps based on pval < 0.05 (for pairwise comparisons)
mat_de_phi_diff_matrix <- mat_de_phi_diff[-c(1)]
mat_de_phi_diff_matrix <- mat_de_phi_diff_matrix/max(mat_de_phi_diff_matrix)
png("../data/Fig6e.png", units="in", width=7, height=7, res=300)
pheatmap(mat_de_phi_diff_matrix,annotation_row = annot_row, annotation_colors = ann_colors, cluster_cols=F,cluster_rows=T,color=colorRampPalette(c("red","white","blue"))(50),show_rownames = F,fontsize = 12)
dev.off()


## DESEQ2 DEgenes vs. GTM Degenes

f <- read.csv(file="../data/Fig5_DEgenes_BRCA_ER_vs_TNBC_together.csv",header=TRUE)
#-log(padj)
png("../data/Fig6d_1.png", units="in", width=7, height=7, res=300)
pheatmap(-log10(f[3]),cluster_rows = F, cluster_cols = F,color=colorRampPalette(c("darkgrey","white"))(50),cellwidth = 20,show_rownames = F,fontsize = 12)
dev.off()
#log2FC
png("../data/Fig6d_2.png", units="in", width=7, height=7, res=300)
pheatmap(f[2],cluster_rows = F, cluster_cols = F,color=colorRampPalette(c("red","white","blue"))(50),cellwidth = 20,show_rownames = F)
dev.off()
#differences in phi
png("../data/Fig6d_3.png", units="in", width=7, height=7, res=300)
pheatmap(as.matrix(f[4:8]),cluster_rows = F, cluster_cols = T,color=colorRampPalette(c("red","white","blue"))(50),show_rownames=F,cellwidth=20,fontsize = 8)
dev.off()