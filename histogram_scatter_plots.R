#setting directory and listing files
###########################################################
setwd("C:/Users/Chester/Desktop/Bonnie/AR81 Excel Sheets")
list.files()
###########################################################

#Visualising the distribution of the expression of a biomarker in a cohort
#Plotting histograms fitted with density curve 
###########################################################

#installing required packages
###########################################################
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("ggthemes")
install.packages("ggrepel")
install.packages("survminer")
###########################################################

#loading the installed packages
###########################################################
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggrepel)
library(ggsci)
library(SummarizedExperiment)
library(readxl)
library(GGally)
library(pander)
library(broom)
library(survminer)
###########################################################

#reading in the dataset and checking the column headers
###########################################################
raw_dat <- read_excel("All TMAs Combined AR81 MBC.xlsx", sheet = "Allred and H-Score Combined")
head(raw_dat)
colnames(raw_dat)
###########################################################

#extracting biomarker expression data 
###########################################################
xpr_dat_hsc <- raw_dat$`Average H-Score`
xpr_dat_alrd <- raw_dat$`Average Allred Score`
###########################################################

#plotting histogram fitted with density curve
###########################################################
pdf("ER_S104 Histogram Average H-Score.pdf", paper = 'a4r')

ggplot(raw_dat, 
       aes(x = xpr_dat_alrd)) + 
  geom_histogram(aes(y = ..density..), 
                 colour = "black", 
                 fill = "white" ) + 
  geom_density(alpha = 0.2, 
               fill = "purple4") + 
  scale_x_continuous(breaks = seq(1, 9, 0.5), 
                     lim = c(1, 9)) + 
  #ylim(0, 0.03) +
  xlab("Average Allred Score") +
  ggtitle("Phosphorylated Androgen Receptor Expression in Male Breast Cancer") +
  theme_few() +
  geom_vline(data = raw_dat, 
             aes(xintercept = median(xpr_dat_alrd)),
             colour = "black",
             linetype = "dashed",
             size = 0.5) +
  geom_text(data = raw_dat,
            aes(label = median(xpr_dat_alrd)),
            x = 7.0,
            y = 0.85,
            size = 5) +
  geom_text(x = 6.5,
            y = 0.85,
            label = "Median = ",
            size = 5)


dev.off()
###########################################################

#checking correlation between replicate average and max biomarker expression using scatter plot and Spearman's correlation
###########################################################
pdf("Average vs Max H-Score Scatter Plot.pdf")

ggplot(raw_dat,
       aes(x = xpr_dat_max,
           y = xpr_dat_avg)) +
  geom_point(shape = 1) +
  xlim(2, 8) +
  ylim(2, 8) +
  stat_cor(method = "spearman",
           label.x = 3,
           label.y = 7) +
  geom_smooth(method = lm,
              colour = "indianred4",
              size = 0.5) +
  xlab("Max Allred Score") +
  ylab("Average Allred Score") +
  theme_few() +
  ggtitle("Average vs Max Allred Scores: Phospho ER Exp in MBC (Serine 294)")

dev.off()
###########################################################

#overlay of multiple biomarker distributions
###########################################################
biomarker <- raw_dat$Biomarker

ggplot(df, 
       aes(x = xpr,
           color = biomarker)) +
  geom_density(size = 1) +
  scale_x_continuous(breaks = seq(1, 9, 0.5), 
                     lim = c(1, 9))+
  xlab("Average Allred Score") +
  theme_few() + 
  scale_color_npg()
###########################################################      
sessionInfo()

