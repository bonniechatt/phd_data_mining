library(readxl)
library(ggplot2)
library(ggpubr)





dat <- read_excel("ESTIMATE_TumorPurity_AllDatasetsCombined.xlsx", sheet = "female_only") 






ggboxplot(dat,
          x = "diff_type_count",
          y = "TumorPurity_Percentage",
          color = "black",
          fill = "diff_type_count",
          xlab = "Number of SSP Models Predicted by genefu",
          ylab = "Tumour Purity",) +
  stat_compare_means(label.x = 0.75, label.y = 150) +
  stat_compare_means(comparisons = diff_types) +
  #geom_jitter(aes(colour = dat$ssp_all)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  ggtitle("Tumour Purity vs Discrepancy in SSP Model Predictions : Female Specific") +
  guides(fill = guide_legend(title = "Number of Models Predicted by genefu"))

diff_types <- list(c(1,2), c(2,3), c(1,3))


dat$ssp_all <- paste(dat$pam50, dat$ssp2006, dat$ssp2003, sep = " + ")

levels(factor(dat$ssp_all))


