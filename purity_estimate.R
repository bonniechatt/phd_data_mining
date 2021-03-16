setwd("C:/Users/Chester/Desktop/Bonnie/data_mining/Callari_data")
list.files()

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos = rforge, dependencies = TRUE)
install.packages("data.table")

library(estimate)
library(data.table)
library(transcripTools)
library(dplyr)

help(package = "estimate")

gset <- fread("S0712.noOutliers-PCATransf.full_list.txt", data.table = F)

gset_final <- gset[, -1]
rownames(gset_final) <- gset$GeneSymbol

eset <- idReplace(gset_final, 
                  format_in = "ensembl_gene_id",
                  format_out = "hgnc_symbol")

eset_final <- eset
eset_final$GeneNames <- rownames(eset)

eset_final <- eset_final[,c(47, 1:46)]

filterCommonGenes(input.f = "cruk_cambridge_mbc.txt", output.f = "Callari_ESTIMATE_common_genes.gct", id = "EntrezID")

estimateScore("Callari_ESTIMATE_common_genes.gct", "Callari_ESTIMATE_scores.gct")

#this only works with Affy
p <- plotPurity(scores = "AlmacDSA_ESTIMATE_scores.gct", samples = "S0712F0014a", platform = "affymetrix")
p

colnames(MBC_Expr[1]) <- "GeneSymbol"


eset_test <- fread("estimate_test_file.txt", data.table = F)
write.table(eset_final, "estimate_test_file.txt", sep = "\t")

eset_test <- eset_test[, -1]
####################################################################################################

MBC_Exp <- system.file("C:/Users/Chester/Desktop/Bonnie/data_mining/GEO_data/GSE104730", 
                       "estimate_test_file.txt", 
                       package = "estimate")

read.table(MBC_Exp)[1:4, 1:4]


####################################################################################################

is.character(eset_test) && length(eset_test) == 1 && nzchar(eset_test)




















