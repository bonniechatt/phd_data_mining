library(mice)
library(readxl)
library(missMDA)
library(sva)
library(factoextra)
library(FactoMineR)
library(sjmisc)

setwd("D:/Missing Gene Imputations")
list.files()

data <- read.delim("combined_dataset_with_gaps_both_sexes.txt", 
                   header = TRUE,
                   sep = "\t")

data <- t(data)
data <- data.frame(data)

data_pam50 <- read.delim("batch_corrected_native_data_both_sexes_pam50.txt",
                         header = TRUE,
                         sep = "\t")

data_pam50 <- t(data_pam50)
data_pam50 <- data.frame(data_pam50)

colnames(data_pam50) <- data_pam50[1,]
data_pam50 <- data_pam50[-1,]

pam50_assign <- data_pam50$pam50_sex
data <- cbind(data, pam50_assign)

batch_combat <- read.delim("combat_batches_both_sexes.txt", 
                           header = TRUE,
                           sep = "\t")
batch_combat <- cbind(batch_combat, pam50_assign)
batch_LumA <- subset(batch_combat, pam50_assign == "Luminal A_Male" | pam50_assign == "Luminal A_Female")
batch_LumB <- subset(batch_combat, pam50_assign == "Luminal B_Male" | pam50_assign == "Luminal B_Female")
batch_HER2 <- subset(batch_combat, pam50_assign == "HER2 Enriched_Male" | pam50_assign == "HER2 Enriched_Female")
batch_Basal <- subset(batch_combat, pam50_assign == "Basal-like_Male" | pam50_assign == "Basal-like_Female")
batch_Normal <- subset(batch_combat, pam50_assign == "Normal-like_Male" | pam50_assign == "Normal-like_Female")


#subsetting data according to pam50 subtypes
data_luminal_a <- subset(data, pam50_assign == "Luminal A_Male" | pam50_assign == "Luminal A_Female")
data_luminal_a$pam50_assign <- NULL 
data_luminal_a <- t(data_luminal_a)
data_luminal_a <- data.frame(data_luminal_a)

data_luminal_b <- subset(data, pam50_assign == "Luminal B_Male" | pam50_assign == "Luminal B_Female")
data_luminal_b$pam50_assign <- NULL
data_luminal_b <- t(data_luminal_b)
data_luminal_b <- data.frame(data_luminal_b)

data_her2_enriched <- subset(data, pam50_assign == "HER2 Enriched_Male" | pam50_assign == "HER2 Enriched_Female")
data_her2_enriched$pam50_assign <- NULL
data_her2_enriched <- t(data_her2_enriched)
data_her2_enriched <- data.frame(data_her2_enriched)

data_basal_like <- subset(data, pam50_assign == "Basal-like_Male" | pam50_assign == "Basal-like_Female")
data_basal_like$pam50_assign <- NULL
data_basal_like <- t(data_basal_like)
data_basal_like <- data.frame(data_basal_like)

data_normal_like <- subset(data, pam50_assign == "Normal-like_Male" | pam50_assign == "Normal-like_Female")
data_normal_like$pam50_assign <- NULL
data_normal_like <- t(data_normal_like)
data_normal_like <- data.frame(data_normal_like)

#imputing data
#nb_lumA <- estim_ncpPCA(data_luminal_a, ncp.min = 0, ncp.max = 5, method.cv = "Kfold", nbsim = 50)
imp_luminal_a <- imputePCA(data_luminal_a,
                           ncp = 2)

imp_luminal_b <- imputePCA(data_luminal_b,
                           ncp = 2)

imp_her2_enriched <- imputePCA(data_her2_enriched,
                           ncp = 2)

imp_basal_like <- imputePCA(data_basal_like,
                           ncp = 2)

imp_normal_like <- imputePCA(data_normal_like,
                           ncp = 2)

#batch correction
imp_lumA_combat <- ComBat(imp_luminal_a$completeObs, batch = batch_LumA$batch)
imp_lumB_combat <- ComBat(imp_luminal_b$completeObs, batch = batch_LumB$batch)
imp_her2_combat <- ComBat(imp_her2_enriched$completeObs, batch = batch_HER2$batch)
imp_basal_combat <- ComBat(imp_basal_like$completeObs, batch = batch_Basal$batch)
imp_normal_combat <- ComBat(imp_normal_like$completeObs, batch = batch_Normal$batch)

#running PCA on imputed dataset
res.pca.LumA <- prcomp(t(imp_lumA_combat))
res.pca.LumB <- prcomp(t(imp_lumB_combat))
res.pca.HER2 <- prcomp(t(imp_her2_combat))
res.pca.Basal <- prcomp(t(imp_basal_combat))
res.pca.Normal <- prcomp(t(imp_normal_combat))

#creating PCA plot for top 20 genes by contribution
#luminal A
pca.lumA <- fviz_pca_ind(res.pca.LumA,
                         legend.title = "Sex",
                         label = "None",
                         geom = "point",
                         pointsize = 3,
                         addEllipses = TRUE,
                         habillage = batch_LumA$pam50_assign,
                         repel = TRUE,
                         title = 'PCA Plot of Imputed Data: Luminal A',
                         palette = c("violetred1", "deepskyblue3")) +
   scale_shape_manual(values = c(19, 19)) 
   
pca.lumA

#luminal B
pca.lumB <- fviz_pca_ind(res.pca.LumB,
                         legend.title = "Sex",
                         label = "None",
                         geom = "point",
                         pointsize = 3,
                         addEllipses = TRUE,
                         habillage = batch_LumB$pam50_assign,
                         repel = TRUE,
                         title = 'PCA Plot of Imputed Data: Luminal B',
                         palette = c("violetred1", "deepskyblue3")) +
   scale_shape_manual(values = c(19, 19)) 

pca.lumB

#her2 enriched
pca.HER2 <- fviz_pca_ind(res.pca.HER2,
                         legend.title = "Sex",
                         label = "None",
                         geom = "point",
                         pointsize = 3,
                         addEllipses = TRUE,
                         habillage = batch_HER2$pam50_assign,
                         repel = TRUE,
                         title = 'PCA Plot of Imputed Data: HER2 Enriched',
                         palette = c("violetred1", "deepskyblue3")) +
   scale_shape_manual(values = c(19, 19)) 

pca.HER2

#basal-like
pca.Basal <- fviz_pca_ind(res.pca.Basal,
                         legend.title = "Sex",
                         label = "None",
                         geom = "point",
                         pointsize = 3,
                         addEllipses = TRUE,
                         habillage = batch_Basal$pam50_assign,
                         repel = TRUE,
                         title = 'PCA Plot of Imputed Data: Basal-like',
                         palette = c("violetred1", "deepskyblue3")) +
   scale_shape_manual(values = c(19, 19)) 

pca.Basal

#normal-like
pca.Normal <- fviz_pca_ind(res.pca.Normal,
                         legend.title = "Sex",
                         label = "None",
                         geom = "point",
                         pointsize = 3,
                         addEllipses = TRUE,
                         habillage = batch_Normal$pam50_assign,
                         repel = TRUE,
                         title = 'PCA Plot of Imputed Data: Normal-like',
                         palette = c("violetred1", "deepskyblue3")) +
   scale_shape_manual(values = c(19, 19)) 

pca.Normal

