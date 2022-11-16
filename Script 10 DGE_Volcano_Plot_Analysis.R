#Script 10

#Original script by Andrew Reinschmidt and adapted by Luis Solano

#Load Libraries#####
library(tidyverse)
library(ggrepel)

#Import DESeq2 generated DEG dataset and save to a variable
deHeLa0Rcnt <-HeLa0RvsHeLaControl_DEG_Results
deHeLa8Rcnt <-HeLa8RvsHeLaControl_DEG_Results
deHeLa8R0R <-HeLa8RvsHeLa0R_DEG_Results

#Create Column of NAs
deHeLa0Rcnt$diffexpressed <- "NO"
deHeLa8Rcnt$diffexpressed <- "NO"
deHeLa8R0R$diffexpressed <- "NO"

#LFC>1.5 and pval<0.05, set as "UP"
deHeLa0Rcnt$diffexpressed[deHeLa0Rcnt$log2FoldChange > 1.5 & deHeLa0Rcnt$padj < 0.05] <- "UP"
deHeLa8Rcnt$diffexpressed[deHeLa8Rcnt$log2FoldChange > 1.5 & deHeLa8Rcnt$padj < 0.05] <- "UP"
deHeLa8R0R$diffexpressed[deHeLa8R0R$log2FoldChange > 1.5 & deHeLa8R0R$padj < 0.05] <- "UP"

#LFC<1.5 and pval<0.05, set as "DOWN"
deHeLa0Rcnt$diffexpressed[deHeLa0Rcnt$log2FoldChange < -1.5 & deHeLa0Rcnt$padj < 0.05] <- "DOWN"
deHeLa8Rcnt$diffexpressed[deHeLa8Rcnt$log2FoldChange < -1.5 & deHeLa8Rcnt$padj < 0.05] <- "DOWN"
deHeLa8R0R$diffexpressed[deHeLa8R0R$log2FoldChange < -1.5 & deHeLa8R0R$padj < 0.05] <- "DOWN"

#Set colors
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

#Create new column "delabel" to de that will contain name of differentially expressed genes
#Gene name not included if diffexpressed = "NO"
deHeLa0Rcnt$delabel <- NA
deHeLa0Rcnt$delabel[deHeLa0Rcnt$diffexpressed != "NO"] <- deHeLa0Rcnt$hgnc_symbol[deHeLa0Rcnt$diffexpressed != "NO"]
deHeLa8Rcnt$delabel <- NA
deHeLa8Rcnt$delabel[deHeLa8Rcnt$diffexpressed != "NO"] <- deHeLa8Rcnt$hgnc_symbol[deHeLa8Rcnt$diffexpressed != "NO"]
deHeLa8R0R$delabel <- NA
deHeLa8R0R$delabel[deHeLa8R0R$diffexpressed != "NO"] <- deHeLa8R0R$hgnc_symbol[deHeLa8R0R$diffexpressed != "NO"]

#Create volcano Plot####
vplotHeLa0RvsHeLaControl<- ggplot(data=deHeLa0Rcnt, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=mycolors) +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

vplotHeLa8RvsHeLaControl<- ggplot(data=deHeLa8Rcnt, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=mycolors) +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

vplotHeLa8RHeLa0R<- ggplot(data=deHeLa8R0R, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=mycolors) +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

#View plots and save to image. Settings below can be changed to export vector based image as well.
vplotHeLa0RvsHeLaControl
ggsave('vplotHeLa0RHeLaControl.png', width =8.50, height =6, device ='png', dpi = 600)

vplotHeLa8RvsHeLaControl
ggsave('vplotHeLa8RHeLaControl.png', width =8.50, height =6, device ='png', dpi = 600)

vplotHeLa8RHeLa0R
ggsave('vplotHeLa8RHeLa0R.png', width =8.50, height =6, device ='png', dpi = 600)
