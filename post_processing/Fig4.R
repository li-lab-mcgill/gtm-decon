## Fig4 - Deconvolution of TCGA-PAAD cohort
setwd("~/Swapna/WORK/McGill/projects/GTM_decon/manuscript/GB_revision/post_processing/")

library(pheatmap)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggpubr)


## Fig4b: Heatmap of TCGA-PAAD with CTS-topics from pancreatic cancer cell types, and de novo bulkRS topics from bulk RNA-seq data
clinical <-read.csv(file="../data/Fig4_TCGA_PAAD_deconvolution.csv")

annot_row <- data.frame("subtype"=clinical$histological_type)
row.names(annot_row) <- row.names(clinical)
annot_row$PNET <- clinical$misclasified
#annot_row$gender <- clinical$gender
#annot_row$stage <- clinical$pathologic_stage
annot_row$days_to_death <- clinical$days_to_death
annot_row$PNET[annot_row$PNET != "Neuroendocrine"] <- "False"
annot_row$PNET[annot_row$PNET == "Neuroendocrine"] <- "True"



#annot_col = data.frame(cell_type = factor(rep(c("Clinical","Pancreas","Immune","Topic"), c(22,14,4,25))))
#annot_col = data.frame(cell_type = factor(rep(c("Clinical","Pancreas","Immune","Topic"), c(22,14,4,10))))
annot_col = data.frame(cell_type = factor(rep(c("Clinical","Pancreas","Topic"), c(22,10,10))))
rownames(annot_col) = colnames(clinical)

#PNET=c(Correct="white",Neuroendocrine="darkred",pseudonormal="yellow",Adjacent_normal_pancreas="yellow",Did_not_arise_from_pancreas="black",Acinar_cell_carcinoma="magenta",Intraductal_papillary_mucinous_neoplasm="cyan",treatment_prior_or_other_malignancy="green")
annot_colors = list(
  subtype=c(pancreas_adenocarcinoma_ductal_type = "#1B9E77", pancreas_adenocarcinoma_other_subtype = "#D95F02", pancreas_colloid_mucinous_non_cystic_carcinoma = "#7570B3", pancreas_undifferentiated_carcinoma = "#E7298A", na = "#FFFFFF"),
  gender=c(female = "pink", male = "blue"),
  #stage=c(stage_i = "#66C2A5",stage_ia = "#FC8D62",stage_ib = "#8DA0CB",stage_iia = "#E78AC3",stage_iib = "#A6D854",stage_iii = "#FFD92F",stage_iv ="#E5C494", na = "#FFFFFF"),
  stage=c(stage_i = "deeppink",stage_ia = "deeppink1",stage_ib = "deeppink2",stage_iia = "deeppink3",stage_iib = "darkviolet",stage_iii = "darkmagenta",stage_iv =" deeppink4", na = "#FFFFFF"),
  T_stage=c(t1 = "#FFFFCC",t2 = "#FED976",t3 = "#FD8D3C",t4 = "#E31A1C",tx = "#FFFFFF", na="#FFFFFF"),
  M_stage=c(m0 = "#FFFFFF",m1 ="#969696",mx = "#FFFFFF"),
  N_stage=c(n0 = "#FFF7FB",n1 = "#D0D1E6", n1b = "#67A9CF",nx = "#FFFFFF", na="#FFFFFF"),
  days_to_death=c("darkgreen","white"),
  cell_type=c(Pancreas = "#0000FF", Topic = "#A9A9A9"),
  lymph_nodes=c("white","green"),
  PNET = c(True = "darkred", False = "white")
)


#Set1 - all topics retained
m <- as.data.frame(select(clinical,23:ncol(clinical)))
#m <- m[-c(10:11,13:14)]
matrix <-t(apply(m,1, function(x) x/sum(x))) #Change each element by value / rowSum

rownames(matrix)<-rownames(clinical)

png("../data/Fig4b.png", units="in", width=10, height=10, res=300)
pheatmap(matrix,annotation_row = annot_row, annotation = annot_col, color=colorRampPalette(c("white", "red"))(50),show_rownames = F,annotation_colors = annot_colors,fontsize_col=12,cellwidth = 15,cex=1)
dev.off()


## Fig4d - Cox-regression for TCGA-PAAD sample cell-type proportions and de novo bulkRS proportions
library("survival")
library("survminer")
library("cluster")
library("factoextra")

#Pancreatic cancer - with immune dataset
d <-read.csv(file="~/Swapna/WORK/McGill/projects/GTM_decon/manuscript/GB_revision/data/Fig4_TCGA_PAAD_deconvolution.csv")

d$days_to_death <- as.numeric(d$days_to_death)
d$days_to_last_followup <- as.numeric(d$days_to_last_followup)
#d <- d[d$days_to_death != "na",]

covariates <- c("years_to_birth", "gender",  "pathologic_stage", "pathology_T_stage", "pathology_N_stage",
                "pathology_M_stage","histological_type","residual_tumor","number_of_lymph_nodes",
                "Acinar","Ductal_type_1","Ductal_type_2","Endocrine","Endothelial","Fibroblast","Stellate","Macrophage",
                "T_cell","B_cell",
                "Topic1","Topic2","Topic3","Topic4","Topic5","Topic6","Topic7",
                "Topic8","Topic9","Topic10")
covariates <- c("Acinar","Ductal_type_1","Ductal_type_2","Endocrine","Endothelial","Fibroblast","Stellate","Macrophage",
                "T_cell","B_cell",
                "Topic1","Topic2","Topic3","Topic4","Topic5","Topic6","Topic7",
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
                   Acinar + Ductal_type_1 + Ductal_type_2 + Endocrine + Endothelial + Fibroblast + Stellate + Macrophage +
                   B_cell + T_cell + Topic1 + Topic2 + Topic3 + Topic4 + Topic5 + Topic6 + Topic7 +
                   Topic8 + Topic9 + Topic10 , data = d)

res.cox <- coxph(Surv(time, vital_status) ~ Acinar + Ductal_type_1 + Ductal_type_2 + Endocrine + Endothelial + Fibroblast + Stellate + Macrophage +
                   B_cell + T_cell + Topic1 + Topic2 + Topic3 + Topic4 + Topic5 + Topic6 + Topic7 +
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


##################### Survival analysis on proportion of Endocrine ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Endocrine), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Endocrine), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Endocrine, 2)
# get cluster means
aggregate(d$Endocrine,by=list(prop_cluster$cluster),FUN=mean)
d$Endocrine_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Endocrine_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Endocrine_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
png("../data/Fig4d.png", units="in", width=10, height=8, res=300)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red", "#2E9FDF"),
           font.x = c(16,"bold"),
           font.y = c(16,"bold"),
           font.tickslab = c(12),
           font.legend = c(14)
)
dev.off()

##################### Survival analysis on proportion of Ductal_type_1 ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Ductal_type_1), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Ductal_type_1), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Ductal_type_1, 2)
# get cluster means
aggregate(d$Ductal_type_1,by=list(prop_cluster$cluster),FUN=mean)
d$Ductal_type_1_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Ductal_type_1_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Ductal_type_1_kmeans, data = d)
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


##################### Survival analysis on proportion of Ductal_type_2 ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Ductal_type_2), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Ductal_type_2), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Ductal_type_2, 2)
# get cluster means
aggregate(d$Ductal_type_2,by=list(prop_cluster$cluster),FUN=mean)
d$Ductal_type_2_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Ductal_type_2_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Ductal_type_2_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
png("../data/KM_ductal_type_2.png", units="in", width=10, height=8, res=300)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red", "#2E9FDF"),
           font.x = c(16,"bold"),
           font.y = c(16,"bold"),
           font.tickslab = c(12),
           font.legend = c(14))
dev.off()

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


##################### Survival analysis on proportion of Stellate ###################
# Determine number of clusters
fviz_nbclust(as.data.frame(d$Stellate), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Stellate), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Stellate, 2)
# get cluster means
aggregate(d$Stellate,by=list(prop_cluster$cluster),FUN=mean)
d$Stellate_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Stellate_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Stellate_kmeans, data = d)
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
fviz_nbclust(as.data.frame(d$Topic5), kmeans, method = "wss")
fviz_nbclust(as.data.frame(d$Topic5), kmeans, method = "silhouette")

#K-means clustering
prop_cluster <- kmeans(d$Topic5, 2)
# get cluster means
aggregate(d$Topic5,by=list(prop_cluster$cluster),FUN=mean)
d$Topic5_kmeans <- prop_cluster$cluster
fit <- survfit(Surv(time,vital_status) ~ Topic5_kmeans, data = d)

#the log rank test for difference in survival 
surv_diff <- survdiff(Surv(time, vital_status) ~ d$Topic5_kmeans, data = d)
surv_diff

#Visualize survival curves
# Change color, linetype by strata, risk.table color by strata
png("../data/KM_Topic5.png", units="in", width=10, height=8, res=300)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red", "#2E9FDF"),
           font.x = c(16,"bold"),
           font.y = c(16,"bold"),
           font.tickslab = c(12),
           font.legend = c(14))
dev.off()