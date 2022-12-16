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


clinical <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts_gtm-decon/COMBINED/PAAD_clinical_notes_pancreas_immune_combined_scaled_i4_k10_UNIQ.csv")

annot_row <- data.frame("subtype"=clinical$histological_type)
row.names(annot_row) <- row.names(clinical)
annot_row$misclassified <- clinical$misclasified
annot_row$gender <- clinical$gender
annot_row$stage <- clinical$pathologic_stage
annot_row$days_to_death <- clinical$days_to_death


#annot_col = data.frame(cell_type = factor(rep(c("Clinical","Pancreas","Immune","Topic"), c(22,14,4,25))))
#annot_col = data.frame(cell_type = factor(rep(c("Clinical","Pancreas","Immune","Topic"), c(22,14,4,10))))
annot_col = data.frame(cell_type = factor(rep(c("Clinical","Pancreas","Topic"), c(22,14,10))))
rownames(annot_col) = colnames(clinical)

annot_colors = list(
  subtype=c(pancreas_adenocarcinoma_ductal_type = "#1B9E77", pancreas_adenocarcinoma_other_subtype = "#D95F02", pancreas_colloid_mucinous_non_cystic_carcinoma = "#7570B3", pancreas_undifferentiated_carcinoma = "#E7298A", na = "#FFFFFF"),
  gender=c(female = "pink", male = "blue"),
  #stage=c(stage_i = "#66C2A5",stage_ia = "#FC8D62",stage_ib = "#8DA0CB",stage_iia = "#E78AC3",stage_iib = "#A6D854",stage_iii = "#FFD92F",stage_iv ="#E5C494", na = "#FFFFFF"),
  stage=c(stage_i = "deeppink",stage_ia = "deeppink1",stage_ib = "deeppink2",stage_iia = "deeppink3",stage_iib = "darkviolet",stage_iii = "darkmagenta",stage_iv =" deeppink4", na = "#FFFFFF"),
  T_stage=c(t1 = "#FFFFCC",t2 = "#FED976",t3 = "#FD8D3C",t4 = "#E31A1C",tx = "#FFFFFF", na="#FFFFFF"),
  M_stage=c(m0 = "#FFFFFF",m1 ="#969696",mx = "#FFFFFF"),
  N_stage=c(n0 = "#FFF7FB",n1 = "#D0D1E6", n1b = "#67A9CF",nx = "#FFFFFF", na="#FFFFFF"),
  days_to_death=c("darkgreen","white"),
  cell_type=c(Pancreas = "#0000FF", Immune = "#00FF00", Topic = "#A9A9A9"),
  lymph_nodes=c("white","green"),
  misclassified=c(Correct="white",Neuroendocrine="darkred",pseudonormal="yellow",Adjacent_normal_pancreas="yellow",Did_not_arise_from_pancreas="black",Acinar_cell_carcinoma="magenta",Intraductal_papillary_mucinous_neoplasm="cyan",treatment_prior_or_other_malignancy="green")
)

#Set1 - all topics retained
m <- as.data.frame(select(clinical,23:ncol(clinical)))
m <- m[-c(10:11,13:14)]
matrix <-t(apply(m,1, function(x) x/sum(x))) #Change each element by value / rowSum

rownames(matrix)<-rownames(clinical)
pheatmap(matrix, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,fontsize = 10)
pheatmap(matrix,annotation_row = annot_row, annotation = annot_col, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 10,fontsize_col=12)

## Order based on fibroblast proportion
ordered_matrix<-matrix[order(-matrix[,16]),]
pheatmap(ordered_matrix,annotation_row = annot_row, annotation = annot_col, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 10,cluster_rows=F)

## Correlation with stage or other params
subset_matrix <- clinical[,c('Acinar','Alpha','Ductal','Fibroblast','Monocyte','T_cell','B_cell')]
matf <- as.data.frame(subset_matrix)
matf$days_to_death <- clinical$days_to_death
matf$stage <- clinical$pathologic_stage

matf$stage <- gsub("ia","i",matf$stage)
matf$stage <- gsub("ib","i",matf$stage)

matf$stage <- factor(matf$stage,levels = c("stage_i","stage_ii","stage_iii","stage_iv"))

matf$T_stage <- clinical$pathology_T_stage
matf$M_stage <- clinical$pathology_M_stage
matf$N_stage <- clinical$pathology_N_stage
matf$T_stage <- factor(matf$T_stage,levels = c("t1","t2","t3","t4"))
matf$N_stage <- factor(matf$N_stage,levels = c("n0","n1","n1b"))
matf$M_stage <- factor(matf$M_stage,levels = c("m0","m1","mx"))

matf$lymph_nodes <- clinical$number_of_lymph_nodes

meltedf<- melt(matf,id.vars = c("stage", "days_to_death","T_stage","N_stage","M_stage","lymph_nodes"))
meltedf <- na.omit(meltedf)

## Comparison with stage
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "stage", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

ggboxplot(meltedf, x="stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

## Comparison with T_stage
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "T_stage", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

compare_means(value ~ T_stage, data = meltedf)
my_comparisons <- list(c("t1","t2"),c("t1","t3"), c("t1","t4"), c("t2","t3"), c("t2", "t4"), c("t3", "t4"))
ggboxplot(meltedf, x="T_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means(comparisons = my_comparisons) 

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

my_comparisons <- list(c("m0","m1"),c("m0","mx"), c("m1","mx"))
ggboxplot(meltedf, x="M_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means(comparisons = my_comparisons) 

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

## bulkRS topics
subset_matrix <- clinical[,c('Topic1','Topic2','Topic3','Topic4','Topic5','Topic6','Topic7','Topic8','Topic9','Topic10')]
matf <- as.data.frame(subset_matrix)
matf$days_to_death <- clinical$days_to_death
matf$stage <- clinical$pathologic_stage

matf$stage <- gsub("ia","i",matf$stage)
matf$stage <- gsub("ib","i",matf$stage)

matf$stage <- factor(matf$stage,levels = c("stage_i","stage_ii","stage_iii","stage_iv"))

matf$T_stage <- clinical$pathology_T_stage
matf$M_stage <- clinical$pathology_M_stage
matf$N_stage <- clinical$pathology_N_stage
matf$T_stage <- factor(matf$T_stage,levels = c("t1","t2","t3","t4"))
matf$N_stage <- factor(matf$N_stage,levels = c("n0","n1","n1b"))
matf$M_stage <- factor(matf$M_stage,levels = c("m0","m1","mx"))

matf$lymph_nodes <- clinical$number_of_lymph_nodes

meltedf<- melt(matf,id.vars = c("stage", "days_to_death","T_stage","N_stage","M_stage","lymph_nodes"))
meltedf <- na.omit(meltedf)

## Comparison with stage
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "stage", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

ggboxplot(meltedf, x="stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

## Comparison with T_stage
my_comparisons <- list(c("t1","t2"),c("t1","t3"), c("t1","t4"), c("t2","t3"), c("t2", "t4"), c("t3", "t4"))
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "T_stage", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

ggboxplot(meltedf, x="T_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

## Comparison with N_stage
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "N_stage", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

ggboxplot(meltedf, x="N_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 


## Comparison with M_stage
ggboxplot(meltedf, x="variable", y="value",  color="variable", add = "jitter", facet.by = "M_stage", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

my_comparisons <- list(c("m0","m1"),c("m0","mx"), c("m1","mx"))
ggboxplot(meltedf, x="M_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means(comparisons = my_comparisons) 

subsetf <- meltedf[meltedf$variable == "Topic5",]
ggboxplot(subsetf, x="M_stage", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 14) + font("ylab",size=14) +
  stat_compare_means() 



## COmparison with number of lymph nodes
ggboxplot(meltedf, x="lymph_nodes", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 

ggboxplot(meltedf, x="days_to_death", y="value",  color="variable", add = "jitter", facet.by = "variable", xlab = "", ylab = "Proportion", ncol=4) + font("xy.text",size = 10) +
  stat_compare_means() 


##Tissue proportion vs. Immune proportion
ti <- as.data.frame(select(clinical,23:36))
im <- as.data.frame(select(clinical,37:40))

tissue_prop <- rowSums(ti)
immune_prop <- rowSums(im)

tissue_immune <- data.frame(tissue_prop,immune_prop)
rownames(tissue_immune) <- rownames(clinical)
rownames(annot_row) <- rownames(clinical)


pheatmap(tissue_immune[order(tissue_immune[,2]),],annotation_row = annot_row, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize = 10,cluster_rows = F,cellwidth = 10)

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

matf$stage <- gsub("ia","i",matf$stage)
matf$stage <- gsub("ib","i",matf$stage)

matf$T_stage <- annot_row$T_stage
matf$M_stage <- annot_row$M_stage
matf$N_stage <- annot_row$N_stage
matf$lymph_nodes <- annot_row$lymph_nodes

meltedf<- melt(matf,id.vars = c("stage", "days_to_death","T_stage","N_stage","M_stage","lymph_nodes"))
meltedf <- na.omit(meltedf)
ggscatter(meltedf, x="days_to_death", y="value",  color="stage", facet.by = "variable")
ggscatter(meltedf, x="days_to_death", y="value",  color="M_stage", facet.by = "variable")
ggscatter(meltedf, x="days_to_death", y="value",  color="T_stage", facet.by = "variable")
ggscatter(meltedf, x="days_to_death", y="value",  color="N_stage", facet.by = "variable")

ggscatter(meltedf, x="lymph_nodes", y="value",  color="stage", facet.by = "variable")
ggscatter(meltedf, x="lymph_nodes", y="value",  color="M_stage", facet.by = "variable")
ggscatter(meltedf, x="lymph_nodes", y="value",  color="T_stage", facet.by = "variable")
ggscatter(meltedf, x="lymph_nodes", y="value",  color="N_stage", facet.by = "variable")

my_comparisons <- list(c("Normal","Basal"),c("Normal","LumA"), c("Normal","LumB"), c("Normal","Her2"), c("Basal", "LumA"), c("Basal", "LumB"), c("Basal", "Her2"), c("LumA","LumB"), c("LumA","Her2"), c("LumB","Her2"))
ggboxplot(meltedf, x="variable", y="value",  color="variable", palette =c("#00AFBB", "#FC4E07"), add = "jitter", facet.by = "stage", xlab = "", ylab = "Combined cell-type proportions per tissue", ncol=4) +
  stat_compare_means(comparisons = my_comparisons)

ggboxplot(meltedf, x="variable", y="value",  color="variable", palette =c("#00AFBB", "#FC4E07"), add = "jitter", facet.by = "T_stage", xlab = "", ylab = "Combined cell-type proportions per tissue") +
  stat_compare_means()

ggboxplot(meltedf, x="variable", y="value",  color="variable", palette =c("#00AFBB", "#FC4E07"), add = "jitter", facet.by = "M_stage", xlab = "", ylab = "Combined cell-type proportions per tissue") +
  stat_compare_means()

ggboxplot(meltedf, x="variable", y="value",  color="variable", palette =c("#00AFBB", "#FC4E07"), add = "jitter", facet.by = "N_stage", xlab = "", ylab = "Combined cell-type proportions per tissue") +
  stat_compare_means()


# b) Kaplan - Meier curves for TCGA PAAD dataset

library("survival")
library("survminer")
library("cluster")
library("factoextra")

#Pancreatic cancer - with immune dataset
d <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/scripts_gtm-decon/COMBINED/PAAD_clinical_notes_pancreas_immune_combined_scaled_i4_k10_UNIQ.csv")

d$days_to_death <- as.numeric(d$days_to_death)
d$days_to_last_followup <- as.numeric(d$days_to_last_followup)
#d <- d[d$days_to_death != "na",]

covariates <- c("years_to_birth", "gender",  "pathologic_stage", "pathology_T_stage", "pathology_N_stage",
                "pathology_M_stage","histological_type","residual_tumor","number_of_lymph_nodes",
                "Acinar","Alpha","Beta","Delta","Ductal","Endothelial","Epsilon","Gamma","Mast",
                "Co.expression","MHC_class_II","PSC","unclassified.endocrine","unclassified",
                "B_cell","Fibroblast","Monocyte","T_cell","Topic1","Topic2","Topic3","Topic4","Topic5","Topic6","Topic7",
                "Topic8","Topic9","Topic10")
covariates <- c("Acinar","Alpha","Beta","Delta","Ductal","Endothelial","Epsilon","Gamma","Mast",
                "Co.expression","MHC_class_II","PSC","unclassified.endocrine","unclassified",
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
res.cox <- coxph(Surv(time, vital_status) ~ years_to_birth + gender + pathologic_stage + pathology_T_stage + pathology_N_stage +
                   pathology_M_stage + histological_type + residual_tumor + number_of_lymph_nodes +
                   Acinar + Alpha + Beta + Delta + Ductal + Endothelial + Epsilon + Gamma + Mast +
                   Co.expression + MHC_class_II + PSC + unclassified.endocrine + unclassified +
                   B_cell + Fibroblast + Monocyte + T_cell + Topic1 + Topic2 + Topic3 + Topic4 + Topic5 + Topic6 + Topic7 +
                   Topic8 + Topic9 + Topic10 , data = d)

res.cox <- coxph(Surv(time, vital_status) ~ Acinar + Alpha + Beta + Delta + Ductal + Endothelial + Epsilon + Gamma + Mast +
                   Co.expression + MHC_class_II + PSC + unclassified.endocrine + unclassified +
                   B_cell + Fibroblast + Monocyte + T_cell + Topic1 + Topic2 + Topic3 + Topic4 + Topic5 + Topic6 + Topic7 +
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


##################### Survival analysis on proportion of Acinar ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Acinar), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Acinar), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Acinar, 2)
# get cluster means
aggregate(d$Acinar,by=list(prop_cluster$cluster),FUN=mean)
d$Acinar_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Acinar_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Acinar_kmeans, data = d)
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


##################### Survival analysis on proportion of Alpha ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Alpha), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Alpha), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Alpha, 2)
# get cluster means
aggregate(d$Alpha,by=list(prop_cluster$cluster),FUN=mean)
d$Alpha_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Alpha_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Alpha_kmeans, data = d)
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


##################### Survival analysis on proportion of Beta ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Beta), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Beta), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Beta, 2)
# get cluster means
aggregate(d$Beta,by=list(prop_cluster$cluster),FUN=mean)
d$Beta_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Beta_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Beta_kmeans, data = d)
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

##################### Survival analysis on proportion of Delta ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Delta), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Delta), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Delta, 2)
# get cluster means
aggregate(d$Delta,by=list(prop_cluster$cluster),FUN=mean)
d$Delta_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Delta_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Delta_kmeans, data = d)
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

##################### Survival analysis on proportion of Ductal ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Ductal), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Ductal), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Ductal, 2)
# get cluster means
aggregate(d$Ductal,by=list(prop_cluster$cluster),FUN=mean)
d$Ductal_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Ductal_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Ductal_kmeans, data = d)
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

##################### Survival analysis on proportion of Endothelial ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Endothelial), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Endothelial), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Endothelial, 2)
# get cluster means
aggregate(d$Endothelial,by=list(prop_cluster$cluster),FUN=mean)
d$Endothelial_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Endothelial_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Endothelial_kmeans, data = d)
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

##################### Survival analysis on proportion of Epsilon ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Epsilon), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Epsilon), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Epsilon, 2)
# get cluster means
aggregate(d$Epsilon,by=list(prop_cluster$cluster),FUN=mean)
d$Epsilon_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Epsilon_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Epsilon_kmeans, data = d)
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

##################### Survival analysis on proportion of Gamma ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Gamma), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Gamma), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Gamma, 2)
# get cluster means
aggregate(d$Gamma,by=list(prop_cluster$cluster),FUN=mean)
d$Gamma_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Gamma_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Gamma_kmeans, data = d)
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

##################### Survival analysis on proportion of PSC ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$PSC), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$PSC), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$PSC, 2)
# get cluster means
aggregate(d$PSC,by=list(prop_cluster$cluster),FUN=mean)
d$PSC_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ PSC_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$PSC_kmeans, data = d)
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
fviz_nbclust(as.data.frame(d$Topic10), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Topic10), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Topic10, 2)
# get cluster means
aggregate(d$Topic10,by=list(prop_cluster$cluster),FUN=mean)
d$Topic10_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Topic10_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Topic10_kmeans, data = d)
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
cox <- coxph(Surv(time, vital_status) ~ Acinar + Alpha + Beta + Delta + Ductal + Endothelial + Epsilon + Gamma + Mast +
               Co.expression + MHC_class_II + PSC + unclassified.endocrine + unclassified +
               B_cell + Fibroblast + Monocyte + T_cell + Topic1 + Topic2 + Topic3 + Topic4 + Topic5 + Topic6 + Topic7 +
               Topic8 + Topic9 + Topic10 , data = d)
summary(cox)
cox_fit <- survfit(cox)
autoplot(cox_fit)

aa_fit <-aareg(Surv(time, vital_status) ~ years_to_birth + gender + pathologic_stage + pathology_T_stage + pathology_N_stage +
                 pathology_M_stage + histological_type + residual_tumor + number_of_lymph_nodes + Alpha + Beta + Delta + Ductal + Endothelial + Epsilon + Gamma + Mast +
                 Co.expression + MHC_class_II + PSC + unclassified.endocrine + unclassified +
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

f <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/TCGA/PAAD/unsup/trainData_JCVB0_nmar_K10_iter99_phi_normalized.csv",header=FALSE)
genes <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/results_gtm-decon/TCGA/PAAD/unsup/genes_bulk",header=FALSE)

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


## Topic5
gseaDat <- df
ranks <- gseaDat$V7
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 

gseaDat <- df
ranks <- gseaDat$V7
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 

gseaDat <- df
ranks <- gseaDat$V7
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(c2_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(c2_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]],
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

## Topic10
gseaDat <- df
ranks <- gseaDat$V12
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(cellMarker_gs, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(cellMarker_gs[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 

gseaDat <- df
ranks <- gseaDat$V12
names(ranks) <- rownames(gseaDat)
fgseaRes <- fgsea(h_gene_sets, ranks, minSize = 15, maxSize = 500, eps = 0.0, scoreType = "pos")
head(fgseaRes[order(pval), ],5)
plotEnrichment(h_gene_sets[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway) + theme(text = element_text(size = 20)) 
