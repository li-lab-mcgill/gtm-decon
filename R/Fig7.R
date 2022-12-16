# Figure : Plot for comparing heatmaps from different reference datasets

library(corrplot)

# a)
f_arti <- read.csv(file="~/WORK/McGill/projects/GTM_decon/scripts/figures/Figure5a_pancreas_cell_proportions_annot.csv",row.names = 1)

#5 topic all genes 
f_seger_1topic <- read.csv(file="~/WORK/McGill/projects/GTM_decon/scripts/figures/Figure5a_PI_Segerstolpe_5topics_pp_artificial_bulkData_trainData_JCVB0_nmar_iter99_metaphe_normalized.csv",header=FALSE)
f_baron_1topic <- read.csv(file="~/WORK/McGill/projects/GTM_decon/scripts/figures/Figure5a_PI_Baron_5topics_pp_artificial_bulkData_trainData_JCVB0_nmar_iter99_metaphe_normalized.csv",header=FALSE)
f_mp_1topic <- read.csv(file="~/WORK/McGill/projects/GTM_decon/scripts/figures/Figure5a_MP_5topics_artificial_bulkData_trainData_JCVB0_nmar_iter99_metaphe_normalized.csv",header=FALSE)

corrplot(as.matrix(f_arti[,c(1:9,12)]),method = 'square', col = COL1("Reds",10), tl.col="black")

colnames(f_seger_1topic)=c("acinar","alpha","beta","delta","ductal","endothelial","epsilon","gamma","mast","co.expression","MHC.classII","PSC","unclassified","unclassified endocrine")
corrplot(as.matrix(f_seger_1topic[,c(1:9,12)]),method = 'square', col = COL1("OrRd",10), tl.col="black")

colnames(f_baron_1topic)=c("acinar","alpha","beta","delta","ductal","endothelial","epsilon","gamma","mast","macrophage","a_stellate","q_stellate","schwann","t_cell")
corrplot(as.matrix(f_baron_1topic[,c(1:9,12)]),method = 'square', col = COL1("Reds",10), tl.col="black")

colnames(f_mp_1topic)=c("alpha","beta","delta","ductal","endothelial","gamma","macrophage","a_stellate","q_stellate","schwann","t_cell","immune_other")
corrplot(as.matrix(f_mp_1topic[,c(1:6,8:9)]),method = 'square', col = COL1("Reds",10), tl.col="black")


# b) Deconvolution of TCGA datasets with different reference profiles

clinical <- read.csv(file="~/WORK/McGill/projects/GTM_decon/scripts/figures/PAAD_clinical_notes.csv")
annot_row <- data.frame("subtype"=clinical$histological_type)
row.names(annot_row) <- row.names(clinical)
annot_row$gender <- clinical$gender
annot_row$stage <- clinical$pathologic_stage
annot_row$days_to_death <- clinical$days_to_death
rgfhgfdxz
annot_row <- data.frame(lapply(annot_row, function(x) {gsub("pancreas-adenocarcinoma ductal type", "pancreas_adenocarcinoma_ductal_type", x)}))
annot_row <- data.frame(lapply(annot_row, function(x) {gsub("pancreas-adenocarcinoma-other subtype", "pancreas_adenocarcinoma_other_subtype", x)}))
annot_row$subtype <- gsub("\\(|\\)", "", annot_row$subtype)
annot_row <- data.frame(lapply(annot_row, function(x) {gsub("pancreas-colloid mucinous non-cystic carcinoma", "pancreas_colloid_mucinous_non_cystic_carcinoma", x)}))
annot_row <- data.frame(lapply(annot_row, function(x) {gsub("pancreas-undifferentiated carcinoma", "pancreas_undifferentiated_carcinoma", x)}))
annot_row[annot_row==""] <- "na"

annot_row <- data.frame(lapply(annot_row, function(x) {gsub("stage ia", "stage i", x)}))
annot_row <- data.frame(lapply(annot_row, function(x) {gsub("stage ib", "stage i", x)}))
annot_row <- data.frame(lapply(annot_row, function(x) {gsub("stage i", "stage_i", x)}))


annot_row$days_to_death <- as.integer(clinical$days_to_death)


# annot_col = data.frame(cell_type = factor(rep(c("Clinical","Pancreas","Immune","Topic"), c(21,14,25,25))))
# annot_col = data.frame(cell_type = factor(rep(c("Clinical","Pancreas", c(21,14,25,25))))
#                        rownames(annot_col) = colnames(clinical)

annot_colors = list(
  subtype=c(pancreas_adenocarcinoma_ductal_type = "#1B9E77", pancreas_adenocarcinoma_other_subtype = "#D95F02", pancreas_colloid_mucinous_non_cystic_carcinoma = "#7570B3", pancreas_undifferentiated_carcinoma = "#E7298A", na = "#FFFFFF"),
  gender=c(female = "#66A61E", male = "#E6AB02"),
  stage=c(stage_i = "#66C2A5",stage_iia = "#E78AC3",stage_iib = "#A6D854",stage_iii = "#FFD92F",stage_iv ="#E5C494", na = "#FFFFFF"),
  T_stage=c(t1 = "#FFFFCC",t2 = "#FED976",t3 = "#FD8D3C",t4 = "#E31A1C",tx = "#FFFFFF", na="#FFFFFF"),
  M_stage=c(m0 = "#FFFFFF",m1 ="#969696",mx = "#FFFFFF"),
  N_stage=c(n0 = "#FFF7FB",n1 = "#D0D1E6", n1b = "#67A9CF",nx = "#FFFFFF", na="#FFFFFF"),
  days_to_death=c("darkgreen","white")
)
#stage=c(stage_i = "#66C2A5",stage_ia = "#FC8D62",stage_ib = "#8DA0CB",stage_iia = "#E78AC3",stage_iib = "#A6D854",stage_iii = "#FFD92F",stage_iv ="#E5C494", na = "#FFFFFF"),
#cell_type=c(Pancreas = "#0000FF", Immune = "#00FF00", Topic = "#A9A9A9")


## Comparison between Human and Mouse reference datasets - for PAAD ##
## Only cell types found in both retained

f <- read.csv(file="~/WORK/McGill/projects/GTM_decon/results_GTM/MP/5topics/pp/TCGA_bulkData_Mouse_trainData_JCVB0_nmar_iter4_metaphe_normalized.csv",header=FALSE)
mat <- as.matrix(f[,c(1:6,8,9)])
colnames(mat)<-c("alpha","beta","delta","ductal","endothelial","gamma","activated_stellate","quiescent_stellate")
matrix <-t(apply(mat,1, function(x) x/sum(x))) #Change each element by value / rowSum
rownames(matrix) <- rownames(f)
mouse_mat <- matrix
hm1<-pheatmap(matrix, annotation_row = annot_row, annotation_colors = annot_colors, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,cluster_cols=F,cluster_rows=T,fontsize = 8, fontsize_col = 11, border_color = NA, legend = TRUE) 

f <- read.csv(file="~/WORK/McGill/projects/GTM_decon/results_GTM/PI_Segerstolpe/all_patients_5topics_0.18each//pp/TCGA_bulkData_trainData_JCVB0_nmar_iter4_metaphe_normalized.csv",header=FALSE)
mat <- as.matrix(f[,c(2:6,8,12)])
colnames(mat)<-c("alpha","beta","delta","ductal","endothelial","gamma","PSC")
matrix <-t(apply(mat,1, function(x) x/sum(x))) #Change each element by value / rowSum
rownames(matrix) <- rownames(f)
reordered_matrix <- matrix[hm1$tree_row$order,] #Re-order on the basis of segerstolpe
pheatmap(reordered_matrix, annotation_row = annot_row, annotation_colors = annot_colors, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,cluster_cols=F,cluster_rows=F,fontsize = 8, fontsize_col = 11, border_color = NA, legend = TRUE) 

