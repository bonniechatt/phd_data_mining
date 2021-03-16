#dependencies
library(mice)
library(VIM)
library(lattice)
library(transcripTools)

#reading in the data frames. Replace names as necessary
dat1 <- read.delim("file1.txt", sep = "\t", header = TRUE)
dat2 <- read.delim("file2.txt", sep = "\t", header = TRUE)
dat3 <- read.delim("file3.txt", sep = "\t", header = TRUE)
dat4 <- read.delim("file4.txt", sep = "\t", header = TRUE)

#setting row names as gene names (column name should be "GeneNames). Replace "dat" with object name
n <- dat
dat <- dat[-1]
rownames(dat) <- n$GeneNames

#if gene ids are in Entrez/Ensembl format, changing them to HGNC symbols. Repeat for each data frames
dat <-idReplace(dat,
                format_in = "entrezgene_id",
                format_out = "hgnc_symbol")

#setting automatic row names and pushing the gene ids into the data set as a separate column identifier 
#this is in preparation of merging the data frames
#assume data frame contains 10,000 genes and 20 patients. Replace as necessary
dat$GeneNames <- rownames(dat)
dat <- dat[,c(21, 1:20)]
rownames(dat) <- 1:10000

#merging data frames
all.male.dat <- merge(dat1, dat2,
                      by = "GeneNames", sort = TRUE, all = TRUE)
all.male.dat <- merge(all.male.dat, dat3,
                      by = "GeneNames", sort = TRUE, all = TRUE)
all.male.dat <- merge(all.male.dat, dat4,
                      by = "GeneNames", sort = TRUE, all = TRUE)


#removing duplicates if any
all.male.dat <- all.male.dat[, !duplicated(colnames(all.male.dat))]

#storing merged data into separate object for analysis and checking summary
mbc.ht <- all.male.dat
summary(mbc.ht)

#setting gene names as row names for the merged dataset
n <- mbc.ht
mbc.ht <- mbc.ht[-1]
rownames(mbc.ht) <- n$GeneNames

#returning table of missing values in each patient
mis.val.num <- md.pattern(mbc.ht)

#visualising missing values
mis.plot <- aggr(mbc.ht, col = c('forestgreen', 'firebrick'),
                  numbers = TRUE, sortVars = TRUE,
                  labels = names(mbc.ht), cex.axis = 0.5,
                  ylab = c("Missing Data", "Pattern"))


#imputing missing values
imp.dump <- mice(mbc.ht, m = 5, maxit = 10, printFlag = TRUE)

#combining imputed values with the merged data frame input
mbc.dat.imp <- mice::complete(imp.dump, "long", inc = TRUE)

#inspecting the distributions of the original and imputed data 
#labels observed data in green and imputed data in red for a given patient
#plotting labels for given patient using stripplot
col <- rep(c("forestgreen", "firebrick")[1 + as.numeric(is.na(imp_dump$data$wz692.tsv))], 6)

stripplot(wz692.tsv~.imp, 
          data = mbc_datimp, 
          jit = TRUE, 
          col = col, 
          xlab = "Number of Imputations")

#not sure if this is the correct way of writing the combined data
write.table(mbc_datimp, "F:/Bonnie/mbc_imputed_data.txt", sep = "\t")

