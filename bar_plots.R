setwd("C:/Users/Chester/Desktop/Bonnie/data_mining")
list.files()

install.packages("readxl")
install.packages("ggplot2")
install.packages("ggsci")
install.packages("circlize")
library(readxl)
library(ggplot2)
library(ggsci)
library(tibble)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggthemes)
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

rm(dat_n)

dat <- read_excel("male_v_female_ssp_compare.xlsx", sheet = "Combined")
head(dat)


dat_n <- within(dat, 
                model <- factor(model, 
                                  levels = names(sort(table(model),
                                                      increasing = TRUE))))

dat_new <- within(dat_n, 
                model <- factor(model, 
                                   levels = names(sort(table(model),
                                                       increasing = TRUE))))
head(dat)

#dat$mol_type <- factor(dat$mol_type, levels = c("Basal-like", "HER2 Enriched", "Normal-like", "Luminal B", "Luminal A"))
dat$model <- factor(dat$model, levels = c("PAM50", "SSP2003", "SSP2006"))



#dat$ssp2006 <- factor(dat$ssp2006, levels = c("Basal-like", "HER2 Enriched", "Normal-like", "Luminal B", "Luminal A"))
#dat$ssp2003 <- factor(dat$ssp2003, levels = c("Basal-like", "HER2 Enriched", "Normal-like", "Luminal B", "Luminal A"))

ggplot(dat, aes(x = model, 
                  fill = mol_type)) +
  geom_bar(position = "fill") +
  #scale_fill_manual(values = c("firebrick", 
                               #"goldenrod", 
                               #"forestgreen")) +
  scale_fill_manual(values = c("firebrick", "palevioletred2", "forestgreen", "deepskyblue2", "midnightblue")) +
  #scale_fill_startrek() +
  ylab("Percentage") +
  xlab("Source Dataset") + 
  #coord_flip() + 
  #scale_y_continuous (labels = scales::percent_format()) +
  #scale_y_continuous(breaks = seq(0, 75, by = 5)) 
  labs(fill = "Molecular Types")

##################################################################################################################
attach(dat)
dat_n <- dat[order(model, mol_type, patient_re_id),]

p <- ggplot(dat, aes(y = sex, x = patient_re_id, fill = molecular_subtype_with_normal)) +
                #fill = mol_type)) +
  #geom_point(cex = 5) +
  geom_tile() +
  #scale_fill_manual(values = c("firebrick", 
  #"goldenrod", 
  #"forestgreen")) +
  scale_fill_manual(values = c("firebrick", "palevioletred2", "midnightblue", "deepskyblue2", "forestgreen")) +
  #scale_fill_startrek() +
  xlab("Patient ID") +
  ylab("Single Sample Predictor Model Type") + 
  theme(axis.text.x = element_text(angle = 90)) +
  labs(fill = "Molecular Types") +
  theme_minimal() + 
  coord_flip()
  #ggtitle("Source Dataset: GSE104730")

p
p + theme(axis.text.x = element_text(angle = 90))

####################################################################################################################
five.mol <- c("Luminal A", "Luminal B", "Normal-like", "HER2 Enriched", "Basal-like")
dat <- dat[dat$mol_type %in% five.mol, ]

dat.scaled <- dat
dat.long <- melt(dat, id = c("patient_re_id", "mol_type", "model"))

heatmap.plot <- ggplot(data = dat.long, aes(x = model, y = patient_re_id)) +
  geom_tile(aes(fill = mol_type)) +
  scale_fill_manual(values = c("firebrick", "palevioletred2", "midnightblue", "deepskyblue2", "forestgreen"))

heatmap.plot


#####################################################################################################################

colors <- structure(1:5, )

cn = colnames(dat)

Heatmap(as.matrix(dat), col = c("firebrick", "palevioletred2", "midnightblue", "deepskyblue2", "forestgreen"),
        cluster_rows = TRUE,
        name = "Molecular Subtypes",
        row_title = "Percentage",
        column_title = "Single Sample Predictor Model (Sex) -- All Datasets Combined",
        row_names_side = "left",
        column_title_side = "bottom",
        #show_column_names = FALSE
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 11),
        rect_gp = gpar(col = "grey", lwd = 0.05))

heatmap(dat)
dat_n <- dat

dat <- dat[, -1]
rownames(dat) <- dat_n$Percentage


