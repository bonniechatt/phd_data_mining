library(mice)
library(readxl)

setwd("D:/")
data <- read_excel("ASR.xlsx", sheet = "Native")

n <- data
data <- data[-1]
rownames(data) <- n$Year

data <- as.data.frame(data)

imp.dump <- mice(data, m = 10, maxit = 10, printFlag = TRUE)

data_imp <- mice::complete(imp.dump, "long", inc = TRUE)

data_imp_broad <- mice::complete(imp.dump, "broad", inc = TRUE)


write.table(data_imp_broad, "D:/ASR_imputed.txt", sep = "\t")
