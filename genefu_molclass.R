setwd("C:/Users/Chester/Desktop/Bonnie")
list.files()

#dependencies
#install.packages("jsonlite")
library(tidyverse)
library(devtools)
library(biomaRt)
library(GEOquery)
library(Biobase)

#devtools::install_github("abc-igmm/transcripTools")
library(transcripTools)
library(genefu)
library(data.table)

#reading in the dataframe
expr_gset <- fread(file = "cruk_cambridge_mbc.csv", data.table = F)

#setting row names
gset_final <- expr_gset[, -1]
row.names(gset_final) <- expr_gset$V1

#replacing NA values with 0
gset_final[is.na(gset_final)] <- 0

#id convert if the format is not entrez id but I'm running it anyway 
#to ensure correctness of format
eset <- idReplace(gset_final, 
                  format_out = "entrezgene_id",
                  format_in = "entrezgene_id")

#source("create_genefu.R")

#getting patient level molecular subtypes
mol.classes <- CreateGenefuPreds(eset)

write.csv(mol.classes, "molclass_cruk_cambridge.csv", quote = F)
