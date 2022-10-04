#Original script by Alison Cher and adapted by Luis Solano
#load libraries
library(ggplot2)
library(reshape2)

#Define boxplot parameters
boxplot.category <- boxplot.category[-c(1), ]
colnames(boxplot.category) <- c('37', '42T0', '42T8')
rownames(boxplot.category) <- c('1','2','3','4','5','6')

boxplot.category_long <- melt(boxplot.category, id=0)
View(boxplot.category_long)
ggplot(boxplot.category_long, aes(x=variable, y=value)) + geom_boxplot()

ggplot(boxplot.formatting, aes(x=Lipid, y=MassSpecValue, color=Treatment)) + geom_boxplot() + geom_jitter(shape = 16, position=position_jitter(0.2))
