#Script 12

#How to annotate and plot PCA
#Need to make a dataframe annotated table that is identical to your 
#numeric vector/data matrix input into prcomp, but w/ text annotations

#You want samples as rows and genes as columns. Last row should be annotations
vsd.df<-as.data.frame(assay(vsd))
head(vsd.df)  

#add a named row to the df above
vsd.df[c('HeatShockStatus'),] <- rbind(c(coldata$HeatShockStatus))
tail(vsd.df)

#transpose to match data format input to PCA
vsdt.df<-t(vsd.df)

#run scaled pca on transpose, autoplot and prepare for real export
vsdt.res<-prcomp(t(assay(vsd)), scale. = TRUE)
HeLa_HeatShockStatus_score<-autoplot(vsdt.res, data = vsdt.df, colour = 'HeatShockStatus')
vsdt.resSummaryObj<-summary(vsdt.res)

ggsave(file="HeLa_HeatShockStatus_score3reps250.tiff",
       plot = HeLa_HeatShockStatus_score,
       device = tiff, 
       units = "in", 
       width=7,
       height = 5,
       dpi = 600)

HeLa_KMC_plot<-autoplot(kmeans(t(assay(vsd)),3), data = t(assay(vsd)))

ggsave(file="HeLa_KMC_plot3reps250.tiff",
       plot = HeLa_KMC_plot,
       device = tiff, 
       units = "in", 
       width=7,
       height = 5,
       dpi = 600)

#Manual kmeans cluster attempt. Start with a transposed vsd normalized df
vsd.dft<-t(vsd.df)
vsd_KMC<-(kmeans(vsd.dft),3))

vsd_KMC$size
vsd_KMC$withinss
summary(vsd_KMC)

vsd.dft$cluster<-factor(vsd_KMC$cluster)
centers<-as.data.frame(vsd_KMC$centers)



#####These aren't working :( #####
#First a standard score plot
pca_HeLa3reps<-ggplot(vsdt.res, aes(x=PC1, y=PC2) ) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",round(100 *vsdt.resSummaryObj$importance[2,1 ]),"% variance")) +
  ylab(paste0("PC2: ",round(100 *vsdt.resSummaryObj$importance[2,2 ]),"% variance")) +
  ggtitle("Score Plot PC1 vs PC2")+
  theme(text = element_text(size = 20))+  
  coord_fixed()

#Then try kmeans clustered and maybe a stat ellipse
pca_HeLa3reps<-ggplot(kmeans(helaVSD_scaledPCA,3), aes(PC1, PC2)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",round(100 *vsdt.resSummaryObj$importance[2,1 ]),"% variance")) +
  ylab(paste0("PC2: ",round(100 *vsdt.resSummaryObj$importance[2,2 ]),"% variance")) +
  ggtitle("Score Plot PC1 vs PC2")+
  theme(text = element_text(size = 20))+  
  coord_fixed()
