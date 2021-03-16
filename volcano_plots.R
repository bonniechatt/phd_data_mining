setwd("C:/Users/Chester/Desktop/Bonnie/Kerri")
list.files()

library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(ggthemes)

genes <- read_excel("probe_results.xlsx")

colnames(genes)
#genes$Significant <- ifelse(genes$log2FoldChange < -(1.5), "Downregulated", 
                            #ifelse(genes$log2FoldChange > 1.5, "Upregulated", "Not Sig"))


genes$Significant <- ifelse(genes$p_value > 0.05, "Not Sig",
                            ifelse(log2(genes$fold_change) < -(0.1), "Downregulated",
                                   ifelse(log2(genes$fold_change) > 0.1, "Upregulated", "Not Sig")))

sig_genes_table <- subset(genes, genes$Significant != "Not Sig")
write.csv(sig_genes_table, "probe_level_sig.csv")

Gene <- genes$gene_id
############################################################################ 
p <- ggplot(genes,
       aes(x = log2(fold_change),
           y = -log10(p_value))) +
  geom_point(alpha = 0.8,
             aes(color = Significant)) +
  geom_vline(
    xintercept = c(-0.1, 0.1),
    col = "gray55",
    linetype = "dotted",
    size = 0.8) +
  scale_color_manual(values = c("seagreen4", "grey", "indianred3")) +
  theme_few(base_size = 12) +
  theme(legend.position = "bottom") +
  geom_text_repel(
    data = subset(genes, -log10(p_value) >= 1.95),
    aes(label = gene_id),
    size = 2.8,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))

p + geom_hline(yintercept = -log10(0.05),
               col = "grey55",
               linetype = "dotted",
               size = 0.8) +
  xlab(expression("Fold Change, Log"[2]*"")) +
  ylab(expression("p-value, Log"[10]*"")) +
  ggtitle("Differential Gene Expression Analysis (Probe Level Analysis):",
          subtitle = "High v. Low Mammographic Density (Median Cut-off)")
#####################################################


data <- read_excel("ERpos_MvDuctalF_Working.xlsx", sheet = "Use_for_ma_plot")

gene_names <- data$gene

data <- data[, -1]
rownames(data) <- make.unique(gene_names)

library(ggpubr)
install.packages("vctrs")

ggmaplot(data,
         main = "Differentially expressed genes in cancer fibroblasts:\nER Positive Male v. Female Invasive Ductal Carcinoma",
         fdr = 0.05,
         fc = 0.5,
         size = 2,
         palette = c("indianred3", "seagreen4", "darkgray"),
         genenames = as.vector(gene_names),
         legend = "top", 
         top = 100,
         font.label = c("bold", 10),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggthemes::theme_few(),
         alpha = 0.5)


#####################################################

exp_data <- read_excel("ERpos_MvDuctalF_Working.xlsx", sheet = "heatmap_desired")

exp_data$log.expression <- log(exp_data$expression)


ggplot(data = exp_data,
       mapping = aes(x = gene,
                     y = subject,
                     fill = log.expression)) + 
  geom_tile() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90,
                                    vjust = 0.5)) +
  scale_fill_viridis()

######################################################

ggplot(data = exp_data,
       mapping = aes(y = gene,
                     x = subject,
                     fill = log.expression)) + 
  geom_tile() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0,
                                   vjust = 0.5)) +
  scale_fill_viridis() +
  ggtitle("Differentially expressed genes in cancer fibroblasts:\nER Positive Male v. Female Invasive Ductal Carcinoma")

######################################################
library(viridis)

table(exp_data$gene)

rm(list = ls(all.names = TRUE))

######################################################

list.files()


library(tidyr)
library(ggplot2)
library(limma)
library(gplots)
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(reshape2)
############################################################################
exp.data <- read_excel("sig_genes_median_absmaxmean.xlsx", sheet = "heatmap")

dat <- exp.data[, 2:409]
rownames(dat) <- exp.data$subject

row.order <- hclust(dist(dat))$order
col.order <- hclust(dist(t(dat)))$order
dat_new <- dat[row.order, col.order]

df_molten_dat <- melt(as.matrix(dat_new))
names(df_molten_dat)[c(1:2)] <- c("Subject", "Gene")

ggplot(data = df_molten_dat,
       aes(x = Subject,
           y = Gene,
           fill = log2(value))) +
  geom_raster()
  



str(exp.data)

exp.long <- pivot_longer(data = exp.data,
                         cols = -c(subject, status),
                         names_to = "gene",
                         values_to = "expression")

exp.long$log.expression <- log(exp.long$expression)

exp.heatmap <- ggplot(data = exp.long,
                      mapping = aes(x = subject,
                                    y = gene,
                                    fill = expression)) +
  geom_tile() +
  scale_fill_gsea() +
  xlab(label = "Patient ID") +
  facet_grid(~ status, 
             switch = "x",
             scales = "free_x",
             space = "free_x") +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5))

exp.heatmap
     







