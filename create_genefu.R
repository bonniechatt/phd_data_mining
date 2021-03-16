CreateGenefuPreds <- function(expressionDF_entrez){
  library(transcripTools)
  library(tidyverse)
  library(genefu)
  
  CreateAnnotationDataframe_entrez <- function(expressionDF){
    df <- data.frame(probe=row.names(expressionDF), EntrezGene.ID=row.names(expressionDF))
    row.names(df) <- df$EntrezGene.ID
    df
  }
  
  CreateAnnotationDataframe_hgnc <- function(expressionDF){
    df <- data.frame(probe=row.names(expressionDF), Gene.Symbol=row.names(expressionDF))
    row.names(df) <- df$Gene.Symbol
    df
  }
  
  an.dfr_entrez <- CreateAnnotationDataframe_entrez(expressionDF = expressionDF_entrez)
  
  #change entrezgene to entrezgene_id
  exp.hgnc <- idReplace(expressionDF_entrez, format_in = "entrezgene_id", "hgnc_symbol")
  an.dfr_hgnc <- CreateAnnotationDataframe_hgnc(exp.hgnc)
  
  
  mammaprint_risk <- gene70(data = t(expressionDF_entrez), annot = an.dfr_entrez, do.mapping = TRUE)$risk %>% as.data.frame()
  rorS_risk <- rorS(data = t(expressionDF_entrez), annot = an.dfr_entrez, do.mapping = TRUE)$risk %>% as.data.frame()
  
  mammaprint_score <- gene70(data = t(expressionDF_entrez), annot = an.dfr_entrez, do.mapping = TRUE)$score %>% as.data.frame()
  rorS_score <- rorS(data = t(expressionDF_entrez), annot = an.dfr_entrez, do.mapping = TRUE)$score %>% as.data.frame()
  
  oncotype_rs <- oncotypedx(data = t(expressionDF_entrez), annot = an.dfr_entrez, do.mapping = TRUE)$risk %>% as.data.frame()
    
  endopredict <- endoPredict(data = t(expressionDF_entrez), annot = an.dfr_entrez, do.mapping = TRUE)$risk %>% as.data.frame()

  
  
#  scmgene <- molecular.subtyping(
#    sbt.model = c("scmgene"), 
#    data = t(expressionDF_entrez), 
#    annot = an.dfr_entrez, do.mapping = TRUE)$subtype %>% 
#    as.data.frame()
  
#  scmod1 <- molecular.subtyping(
#    sbt.model = c("scmod1"), 
#    data = t(expressionDF_entrez), 
#    annot = an.dfr_entrez, do.mapping = TRUE)$subtype %>% 
#    as.data.frame()
  
#  scmod2 <- molecular.subtyping(
#    sbt.model = c("scmod2"), 
#    data = t(expressionDF_entrez), 
#    annot = an.dfr_entrez, do.mapping = TRUE)$subtype %>% 
#    as.data.frame()
  
  pam50 <- molecular.subtyping(
    sbt.model = c("pam50"), 
    data = t(expressionDF_entrez), 
    annot = an.dfr_entrez, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  ssp2006 <- molecular.subtyping(
    sbt.model = c("ssp2006"), 
    data = t(expressionDF_entrez), 
    annot = an.dfr_entrez, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  ssp_2003 <- molecular.subtyping(
    sbt.model = c("ssp2003"), 
    data = t(expressionDF_entrez), 
    annot = an.dfr_entrez, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  intClust <- molecular.subtyping(
    sbt.model = c("intClust"), 
    data = t(exp.hgnc), 
    annot = an.dfr_hgnc, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
 # AIMS <- molecular.subtyping(
    #sbt.model = c("AIMS"), 
    #data = t(expressionDF_entrez), 
    #annot = an.dfr_entrez, do.mapping = TRUE)$subtype %>% 
    #as.data.frame()
  
#  subtypes <- cbind(scmgene, scmod1, scmod2, pam50, ssp2006, ssp_2003, intClust, mammaprint_risk, mammaprint_score, rorS_risk, rorS_score, oncotype_rs, endopredict)
#  colnames(subtypes) <- c("scmgene", "smcod1", "scmod2", "pam50", "ssp2006", "ssp2003", "IC10", "MammaPrint_risk", "MammaPrint_score", "rorS_risk", "rorS_score","oncotype_rs", "endopredict")
  subtypes <- cbind(pam50, ssp2006, ssp_2003, intClust, mammaprint_risk, mammaprint_score, rorS_risk, rorS_score, oncotype_rs, endopredict)
  colnames(subtypes) <- c("pam50", "ssp2006", "ssp2003", "IC10", "MammaPrint_risk", "MammaPrint_score", "rorS_risk", "rorS_score","oncotype_rs", "endopredict")
  subtypes
}

