setwd("C:/Users/Chester/Desktop/Bonnie/data_mining")

BiocManager::install("GEOquery")
BiocManager::install("Biobase")

library(GEOquery)
library(Biobase)

gset <- getGEO("GSE1477", GSEMatrix = TRUE, getGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL1198", attr(gset, "names")) else idx <- 1

gset <- gset[[idx]]

rm(list = ls(all.names = TRUE))

pheno <- pData(gset)

write.table(pheno, file = "GSE1477_metadata.txt", sep = "\t")
