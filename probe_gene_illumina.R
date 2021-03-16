rm(list = ls(all.names = TRUE))

listAttributes(ensembl)


ensembl<-useMart("ensembl", dataset="hsapiens_gene_ensembl")

genes <- rownames(expr_gset)

df <- getBM(attributes = c("hgnc_symbol", "illumina_humanht_12_v3"), 
                    filters = c("illumina_humanht_12_v3"),
                    values = genes, 
                    mart = ensembl)

genes <- merge(x = as.data.frame(genes),
               y = df, 
               by.y="illumina_humanht_12_v3", 
               all.x=T,
               by.x="genes")

genes <- as.matrix(genes, row.names = TRUE)

#need to run this command line twice. The second time returns the correct results without NA values for reasons I don't completely understand
###########################################################
df.final <- getBM(attributes = c("hgnc_symbol", "illumina_humanht_12_v3"), 
                          filters = c("illumina_humanht_12_v3"),
                          values = genes, 
                          mart = ensembl)

head(df.final)


m <- expr_gset
head(m)
expr_gset_nadrop <- m[apply(m, 1, Compose(is.finite, all)),]


gene_list <- df.final$hgnc_symbol
probe_list <- df.final$illumina_humanht_12_v3

keep = !duplicated(probe_list)
allProbes = probe_list[keep]
allGenes = gene_list[keep]

PR = allProbes
GE = allGenes

head(PR)
head(GE)

collapsed_expr_maxmean <- collapseRows(expr_gset_nadrop, GE, PR, method = "MaxMean")

eset <- collapsed_expr_maxmean$datETcollapsed

