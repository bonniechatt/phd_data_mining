library(readxl)
install.packages("wesanderson")
library(wesanderson)
library(ggsci)

getwd()
setwd("C:/Users/Chester/Desktop/Bonnie/ER")
list.files()

df <- read_excel("For Density Plot Overlay.xlsx")
head(df)

colnames(df)

xpr <- df$`Average Allred Score`
ypr <- df$`Average H-Score`
biomarker <- df$Biomarker


ggplot(df, 
       aes(x = xpr,
           color = biomarker)) +
  geom_density(size = 1) +
  scale_x_continuous(breaks = seq(1, 9, 0.5), 
                     lim = c(1, 9))+
  xlab("Average Allred Score") +
  theme_few() + 
  scale_color_npg()


