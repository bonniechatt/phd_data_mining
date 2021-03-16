#setting directory and listing files
###########################################################
setwd("C:/Users/Chester/Desktop/Bonnie/ER/ER_S118 Excel Sheets")
list.files()
###########################################################

#Kaplan Meier plots created using "survminer" 
#Based on both median and data-driven dichotomisation
###########################################################

#installing required non-bioconductor packages
###########################################################
install.packages("survival")
install.packages("magrittr")
install.packages("ggplot2")
install.packages("GGally")
install.packages("ggthemes")
install.packages("cowplot")
install.packages("knitr")
install.packages("broom")
install.packages("pander")
install.packages("readr")
install.packages("rmarkdown")
install.packages("survsim")
install.packages("testthat")
install.packages("desiR")
install.packages("viridis")
install.packages("readxl")
install.packages("cgdsr")
install.packages("convert")
###########################################################

#installing required bioconductor packages
###########################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("survcomp")
BiocManager::install("Biobase")
BiocManager::install("convert")
###########################################################

#installing "survivALL" package from local drive as its absent from CRAN
###########################################################
install.packages("C:/Users/Chester/Desktop/Bonnie/survivALL/survivALL_0.9.3.tar.gz", repos = NULL, type = "source")
###########################################################

#loading the installed packages 
###########################################################
library(survivALL)
library(survival)
library(magrittr)
library(ggplot2)
library(GGally)
library(ggthemes)
library(cowplot)
library(knitr)
library(broom)
library(pander)
library(readr)
library(rmarkdown)
library(survsim)
library(testthat)
library(desiR)
library(viridis)
library(readxl)
library(survcomp)
library(Biobase)
library(cgdsr)
library(convert)
###########################################################

#reading in the dataset and checking the column headers
###########################################################
disc <- read_excel("MBC IHC Matt - Working Table.xlsx", sheet = "OS Data")
head(disc)
colnames(disc)
###########################################################

#creating a median driven binary classifier 
###########################################################
xpr <- disc$`Average H-Score`
med_classifier <- ifelse(xpr >= median(xpr), "high", "low") 
###########################################################

#creating survival object based on the median dichotomisation
###########################################################
srv_obj <- survival::Surv(disc$`OS Time`, disc$OSCALL)
med_fit <- survival::survfit(srv_obj ~ med_classifier)
###########################################################

#plotting KM curve based on median dichotomisation
###########################################################
ggsurvplot(med_fit,
           data = disc,
           #palette = c("indianred3", "darkseagreen4"),
           palette = "ucscgb",
           legend.title = "Expression Category",
           legend.labs = c("High", "Low"),
           size = 0.8,
           pval = TRUE,
           pval.coord = c(250, 0.95),
           pval.method = TRUE,
           pval.method.coord = c(250, 0.9),
           risk.table = TRUE,
           risk.table.y.text.col = TRUE,
           ggtheme = theme_few(),
           xlim = c(0, 300),
           break.time.by = 60,
           title = "Patient stratification by median dichotomisation of \nphosphoAR Avg Allred Score (OS)",
           xlab = "Months",
           font.x = 12,
           font.y = 12,
           font.tickslab = 10,
           font.title = 14)
###########################################################

#checking univariate Cox PH results of the median-driven analysis
###########################################################
summary(survival::coxph(srv_obj ~ med_classifier))
###########################################################

#determining the best point-of-separation using the "survivALL()" function
#creating a binary classifier based on the best point-of-separation determined
###########################################################
srvall <- survivALL(measure = xpr,
                    srv = disc,
                    time = "OS Time",
                    event = "OSCALL",
                    measure_name = "phosphoER")

srvall

data_classifier <- ifelse(srvall$clsf == 0, "low", "high")
#data_classifier <- ifelse(srvall$index >= 9, "high", "low")
###########################################################

#visualising the best point-of-separation using the "plotALL()" function
###########################################################
median_index_val <- nrow(srvall)/2

optimal_index_val <- which(diff(srvall$clsf) == 1)
#optimal_index_val <- which(srvall$index == 72)

annotation_dfr <- data.frame(index = c(median_index_val, optimal_index_val),
                             y = c(-0.2, -1.35),
                             annot = c("Median cut-off\n|", "|\nData-driven cut-off"))

plotALL(measure = xpr,
        srv = disc,
        time = "OS Time",
        event = "OSCALL",
        measure_name = "phosphoAR") +
  geom_text(aes(x = index, 
                y = y, 
                label = annot), 
            data = annotation_dfr) +
  labs(title = "Best point-of-separation for data-driven dichotomisation \nof phosphoAR (Average Allred OS)",
       subtitle = "survivALL plot output: Hazard ratios and significance are combined to calculate desirability (colour scale)\nMedian and data-driven best point-of-separation are labelled") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlim(0, 138) +
  ylim(-2.5, 2.5)
  
###########################################################

#creating a survival object based on the data-driven binary classifier
###########################################################
srv_obj <- survival::Surv(as.numeric(srvall$event_time), srvall$event)
data_fit <- survival::survfit(srv_obj ~ data_classifier)
###########################################################

#plotting KM curve based on data driven dichotomisation
###########################################################
pdf()

ggsurvplot(data_fit,
           data = disc,
           palette = "ucscgb",
           legend.title = "Expression Category",
           legend.labs = c("High", "Low"),
           size = 0.8,
           pval = TRUE,
           pval.coord = c(250, 1.00),
           pval.method = TRUE,
           pval.method.coord = c(250, 0.95),
           risk.table = TRUE,
           risk.table.y.text.col = TRUE,
           ggtheme = theme_few(),
           xlim = c(0, 300),
           break.time.by = 60,
           title = "Patient stratification by data-driven dichotomisation of \nphosphoER_S294 avg Allred Score (OS)",
           xlab = "Months",
           font.x = 12,
           font.y = 12,
           font.tickslab = 10,
           font.title = 14)

dev.off()
###########################################################

#checking univariate Cox PH results of the data-driven analysis
###########################################################
summary(survival::coxph(srv_obj ~ data_classifier))
###########################################################

sessionInfo()
