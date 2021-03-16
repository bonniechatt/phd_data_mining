BiocManager::install("ArrayExpress")

install.packages("pdfInfoBuilder")
library(ArrayExpress)
library(arrayQualityMetrics)
library(geneplotter)
library(vsn)
library(affy)
library(limma)
library(marray)

setwd("C:/Users/Chester/Desktop/Bonnie/Data Mining")

rm(list = ls(all.names = TRUE))

eset <- try(ArrayExpress("E-TABM-810"))
 

