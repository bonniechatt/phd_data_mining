setwd("C:/Users/Chester/Desktop/Bonnie")

list.files()

library(tidyverse)
library(devtools)
library(biomaRt)
library(GEOquery)
library(Biobase)
install.packages("fansi")

devtools::install_github("abc-igmm/transcripTools")
BiocManager::install("genefu")

library(transcripTools)
library(genefu)
library(readxl)
library(data.table)

expr_gset <- fread("cruk_cambridge_mbc.csv", data.table = F)
#expr_gset <- read_excel("MBC_data_collapsed_to_symbols.gct.xlsx", sheet = "to_analyse")
df <- read.delim("S0712.noOutliers-PCATransf.full_list.txt", sep = "\t", header = TRUE)

df$ensembl_gene_edit <- substr(df$Ensembl.Gene.ID, 1, 15)

class(expr_gset)

##################################################

m <- expr_gset
head(m)
expr_gset_nadrop <- m[apply(m, 1, Compose(is.finite, all)),]


gene_list <- df$ensembl_gene_edit
probe_list <- df$Probe.Set.ID

keep = !duplicated(probe_list)
allProbes = probe_list[keep]
allGenes = gene_list[keep]

PR = allProbes
GE = allGenes

head(PR)
head(GE)

#gset <- getGEO("GSE31259", GSEMatrix = TRUE, getGPL = FALSE)

#if (length(gset) > 1) idx <- grep("GPL6947", attr(gset, "names")) else idx <- 1

#gset <- gset[[idx]]

#expr_gset <- exprs(gset)
###################################################
#gene.names <- expr_gset[, 1]
#expr_gset_final <- expr_gset[, -1]
#expr_gset <- expr_gset[, -1]
#row.names(expr_gset_final) <- expr_gset$entrez_id 
###################################################

eset <- idReplace(gset_final, 
               format_in = "entrezgene_id",
               format_out = "entrezgene_id")

source("create_genefu.R")

mol.classes <- CreateGenefuPreds(eset)

write.csv(mol.classes, "molclass_pancancer_2018.csv", quote = F)

gset_final <- gset
gset_final[is.na(gset_final)] <- 1













