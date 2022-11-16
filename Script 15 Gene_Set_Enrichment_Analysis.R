#script 15#
#This script reproduces the functional enrichment analysis process for all HeLa samples with clusterProfiler and gProfiler's GSEA functionalities

#Load libraries#####
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(DT)
library(gprofiler2)
library(enrichplot)

#Create local internal collection objects with msigdbr. #####
#From specified human specific collections, only gene set names and comprising genes are saved

hs_gsea_H <- msigdbr(species = "Homo sapiens", 
                     category = "H") %>% 
  dplyr::select(gs_name, human_ensembl_gene)
hs_gsea_C2 <- msigdbr(species = "Homo sapiens", 
                      category = "C2") %>% 
  dplyr::select(gs_name, human_ensembl_gene)
hs_gsea_C4 <- msigdbr(species = "Homo sapiens", 
                      category = "C4") %>% 
  dplyr::select(gs_name, human_ensembl_gene)
hs_gsea_C5 <- msigdbr(species = "Homo sapiens", 
                      category = "C5") %>% 
  dplyr::select(gs_name, human_ensembl_gene)
hs_gsea_C6 <- msigdbr(species = "Homo sapiens", 
                      category = "C6") %>% 
  dplyr::select(gs_name, human_ensembl_gene)

#Create GSEA input data structures #####

#Create pre-GSEA objects (0RvsControl)
HeLa0RvsHeLaControl.inputstx <-as.data.frame(HeLa0RvsHeLaControl@listData$log2FoldChange, row.names = HeLa0RvsHeLaControl@rownames) %>% tibble::rownames_to_column("GeneID") %>% `colnames<-`(c("GeneID","log2FoldChange"))
HeLa0RvsHeLaControl.gsea<-HeLa0RvsHeLaControl.inputstx$log2FoldChange
names(HeLa0RvsHeLaControl.gsea)<- HeLa0RvsHeLaControl.inputstx$GeneID
HeLa0RvsHeLaControl.gsea <- sort(HeLa0RvsHeLaControl.gsea, decreasing = TRUE)

#Create pre-GSEA objects (8RvsControl)
HeLa8RvsHeLaControl.inputstx <-as.data.frame(HeLa8RvsHeLaControl@listData$log2FoldChange, row.names = HeLa8RvsHeLaControl@rownames) %>% tibble::rownames_to_column("GeneID") %>% `colnames<-`(c("GeneID","log2FoldChange"))
HeLa8RvsHeLaControl.gsea<-HeLa8RvsHeLaControl.inputstx$log2FoldChange
names(HeLa8RvsHeLaControl.gsea)<- HeLa8RvsHeLaControl.inputstx$GeneID
HeLa8RvsHeLaControl.gsea <- sort(HeLa8RvsHeLaControl.gsea, decreasing = TRUE)

#Create pre-GSEA objects (8Rvs0R)
HeLa8RvsHeLa0R.inputstx <-as.data.frame(HeLa8RvsHeLa0R@listData$log2FoldChange, row.names = HeLa8RvsHeLa0R@rownames) %>% tibble::rownames_to_column("GeneID") %>% `colnames<-`(c("GeneID","log2FoldChange"))
HeLa8RvsHeLa0R.gsea<-HeLa8RvsHeLa0R.inputstx$log2FoldChange
names(HeLa8RvsHeLa0R.gsea)<- HeLa8RvsHeLa0R.inputstx$GeneID
HeLa8RvsHeLa0R.gsea <- sort(HeLa8RvsHeLa0R.gsea, decreasing = TRUE)

#Run GSEA#####

#Run and store GSEA result against collections of interest (0RvsControl)
HeLa0RvsHeLaControl.gseaH.res <- GSEA(HeLa0RvsHeLaControl.gsea, TERM2GENE=hs_gsea_H, verbose=TRUE)
HeLa0RvsHeLaControl.gseaC2.res <- GSEA(HeLa0RvsHeLaControl.gsea, TERM2GENE=hs_gsea_C2, verbose=TRUE)
HeLa0RvsHeLaControl.gseaC4.res <- GSEA(HeLa0RvsHeLaControl.gsea, TERM2GENE=hs_gsea_C4, verbose=TRUE)
HeLa0RvsHeLaControl.gseaC5.res <- GSEA(HeLa0RvsHeLaControl.gsea, TERM2GENE=hs_gsea_C5, verbose=TRUE)
HeLa0RvsHeLaControl.gseaC6.res <- GSEA(HeLa0RvsHeLaControl.gsea, TERM2GENE=hs_gsea_C6, verbose=TRUE)

#Run and store GSEA result against collections of interest (8RvsControl)
HeLa8RvsHeLaControl.gseaH.res <- GSEA(HeLa8RvsHeLaControl.gsea, TERM2GENE=hs_gsea_H, verbose=TRUE)
HeLa8RvsHeLaControl.gseaC2.res <- GSEA(HeLa8RvsHeLaControl.gsea, TERM2GENE=hs_gsea_C2, verbose=TRUE)
HeLa8RvsHeLaControl.gseaC4.res <- GSEA(HeLa8RvsHeLaControl.gsea, TERM2GENE=hs_gsea_C4, verbose=TRUE)
HeLa8RvsHeLaControl.gseaC5.res <- GSEA(HeLa8RvsHeLaControl.gsea, TERM2GENE=hs_gsea_C5, verbose=TRUE)
HeLa8RvsHeLaControl.gseaC6.res <- GSEA(HeLa8RvsHeLaControl.gsea, TERM2GENE=hs_gsea_C6, verbose=TRUE)

#Run and store GSEA result against collections of interest (8Rvs0R)
HeLa8RvsHeLa0R.gseaH.res <- GSEA(HeLa8RvsHeLa0R.gsea, TERM2GENE=hs_gsea_H, verbose=TRUE)
HeLa8RvsHeLa0R.gseaC2.res <- GSEA(HeLa8RvsHeLa0R.gsea, TERM2GENE=hs_gsea_C2, verbose=TRUE)
HeLa8RvsHeLa0R.gseaC4.res <- GSEA(HeLa8RvsHeLa0R.gsea, TERM2GENE=hs_gsea_C4, verbose=TRUE)
HeLa8RvsHeLa0R.gseaC5.res <- GSEA(HeLa8RvsHeLa0R.gsea, TERM2GENE=hs_gsea_C5, verbose=TRUE)
HeLa8RvsHeLa0R.gseaC6.res <- GSEA(HeLa8RvsHeLa0R.gsea, TERM2GENE=hs_gsea_C6, verbose=TRUE)

#####Uncomment lines below to export collection specific GSEA results to CSV files

#write.csv(HeLa0RvsHeLaControl.gseaH.res, file = "HeLa0RvsHeLaControl.gseaH.res_NovogeneRawCounts.csv")
#write.csv(HeLa0RvsHeLaControl.gseaC2.res, file = "HeLa0RvsHeLaControl.gseaC2.res_NovogeneRawCounts.csv")
#write.csv(HeLa0RvsHeLaControl.gseaC4.res, file = "HeLa0RvsHeLaControl.gseaC4.res_NovogeneRawCounts.csv")
#write.csv(HeLa0RvsHeLaControl.gseaC5.res, file = "HeLa0RvsHeLaControl.gseaC5.res_NovogeneRawCounts.csv")
#write.csv(HeLa0RvsHeLaControl.gseaC6.res, file = "HeLa0RvsHeLaControl.gseaC6.res_NovogeneRawCounts.csv")

#write.csv(HeLa8RvsHeLaControl.gseaH.res, file = "HeLa8RvsHeLaControl.gseaH.res_NovogeneRawCounts.csv")
#write.csv(HeLa8RvsHeLaControl.gseaC2.res, file = "HeLa8RvsHeLaControl.gseaC2.res_NovogeneRawCounts.csv")
#write.csv(HeLa8RvsHeLaControl.gseaC4.res, file = "HeLa8RvsHeLaControl.gseaC4.res_NovogeneRawCounts.csv")
#write.csv(HeLa8RvsHeLaControl.gseaC5.res, file = "HeLa8RvsHeLaControl.gseaC5.res_NovogeneRawCounts.csv")
#write.csv(HeLa8RvsHeLaControl.gseaC6.res, file = "HeLa8RvsHeLaControl.gseaC6.res_NovogeneRawCounts.csv")

#write.csv(HeLa8RvsHeLa0R.gseaH.res, file = "HeLa8RvsHeLa0R.gseaH.res_NovogeneRawCounts.csv")
#write.csv(HeLa8RvsHeLa0R.gseaC2.res, file = "HeLa8RvsHeLa0R.gseaC2.res_NovogeneRawCounts.csv")
#write.csv(HeLa8RvsHeLa0R.gseaC4.res, file = "HeLa8RvsHeLa0R.gseaC4.res_NovogeneRawCounts.csv")
#write.csv(HeLa8RvsHeLa0R.gseaC5.res, file = "HeLa8RvsHeLa0R.gseaC5.res_NovogeneRawCounts.csv")
#write.csv(HeLa8RvsHeLa0R.gseaC6.res, file = "HeLa8RvsHeLa0R.gseaC6.res_NovogeneRawCounts.csv")

#####Generate Dot plot input tibbles 
#Generate Dot plot input tibbles (0RvsControl)
HeLa0RvsHeLaControl.GSEAHresult<-as_tibble(HeLa0RvsHeLaControl.gseaH.res@result)
HeLa0RvsHeLaControl.GSEAHresult.tib <- HeLa0RvsHeLaControl.GSEAHresult %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA0R",
    NES < 0 ~ "HELAcontrol"))

HeLa0RvsHeLaControl.GSEAC2result<-as_tibble(HeLa0RvsHeLaControl.gseaC2.res@result)
HeLa0RvsHeLaControl.GSEAC2result.tib <- HeLa0RvsHeLaControl.GSEAC2result %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA0R",
    NES < 0 ~ "HELAcontrol"))

HeLa0RvsHeLaControl.GSEAC4result<-as_tibble(HeLa0RvsHeLaControl.gseaC4.res@result)
HeLa0RvsHeLaControl.GSEAC4result.tib <- HeLa0RvsHeLaControl.GSEAC4result %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA0R",
    NES < 0 ~ "HELAcontrol"))

HeLa0RvsHeLaControl.GSEAC5result<-as_tibble(HeLa0RvsHeLaControl.gseaC5.res@result)
HeLa0RvsHeLaControl.GSEAC5result.tib <- HeLa0RvsHeLaControl.GSEAC5result %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA0R",
    NES < 0 ~ "HELAcontrol"))

HeLa0RvsHeLaControl.GSEAC6result<-as_tibble(HeLa0RvsHeLaControl.gseaC6.res@result)
HeLa0RvsHeLaControl.GSEAC6result.tib <- HeLa0RvsHeLaControl.GSEAC6result %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA0R",
    NES < 0 ~ "HELAcontrol"))

#Generate Dot plot input tibbles (8RvsControl)
HeLa8RvsHeLaControl.GSEAHresult<-as_tibble(HeLa8RvsHeLaControl.gseaH.res@result)
HeLa8RvsHeLaControl.GSEAHresult.tib <- HeLa8RvsHeLaControl.GSEAHresult %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA8R",
    NES < 0 ~ "HELAcontrol"))

HeLa8RvsHeLaControl.GSEAC2result<-as_tibble(HeLa8RvsHeLaControl.gseaC2.res@result)
HeLa8RvsHeLaControl.GSEAC2result.tib <- HeLa8RvsHeLaControl.GSEAC2result %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA8R",
    NES < 0 ~ "HELAcontrol"))

HeLa8RvsHeLaControl.GSEAC4result<-as_tibble(HeLa8RvsHeLaControl.gseaC4.res@result)
HeLa8RvsHeLaControl.GSEAC4result.tib <- HeLa8RvsHeLaControl.GSEAC4result %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA8R",
    NES < 0 ~ "HELAcontrol"))

HeLa8RvsHeLaControl.GSEAC5result<-as_tibble(HeLa8RvsHeLaControl.gseaC5.res@result)
HeLa8RvsHeLaControl.GSEAC5result.tib <- HeLa8RvsHeLaControl.GSEAC5result %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA8R",
    NES < 0 ~ "HELAcontrol"))

HeLa8RvsHeLaControl.GSEAC6result<-as_tibble(HeLa8RvsHeLaControl.gseaC6.res@result)
HeLa8RvsHeLaControl.GSEAC6result.tib <- HeLa8RvsHeLaControl.GSEAC6result %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA8R",
    NES < 0 ~ "HELAcontrol"))

#Generate Dot plot input tibbles (8Rvs0R)
HeLa8RvsHeLa0R.GSEAHresult<-as_tibble(HeLa8RvsHeLa0R.gseaH.res@result)
HeLa8RvsHeLa0R.GSEAHresult.tib <- HeLa8RvsHeLa0R.GSEAHresult %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA8R",
    NES < 0 ~ "HELA0R"))

HeLa8RvsHeLa0R.GSEAC2result<-as_tibble(HeLa8RvsHeLa0R.gseaC2.res@result)
HeLa8RvsHeLa0R.GSEAC2result.tib <- HeLa8RvsHeLa0R.GSEAC2result %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA8R",
    NES < 0 ~ "HELA0R"))

HeLa8RvsHeLa0R.GSEAC4result<-as_tibble(HeLa8RvsHeLa0R.gseaC4.res@result)
HeLa8RvsHeLa0R.GSEAC4result.tib <- HeLa8RvsHeLa0R.GSEAC4result %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA8R",
    NES < 0 ~ "HELA0R"))

HeLa8RvsHeLa0R.GSEAC5result<-as_tibble(HeLa8RvsHeLa0R.gseaC5.res@result)
HeLa8RvsHeLa0R.GSEAC5result.tib <- HeLa8RvsHeLa0R.GSEAC5result %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA8R",
    NES < 0 ~ "HELA0R"))

HeLa8RvsHeLa0R.GSEAC6result<-as_tibble(HeLa8RvsHeLa0R.gseaC6.res@result)
HeLa8RvsHeLa0R.GSEAC6result.tib <- HeLa8RvsHeLa0R.GSEAC6result %>%
  mutate(Phenotype = case_when(
    NES > 0 ~ "HELA8R",
    NES < 0 ~ "HELA0R"))

#Create Dotplot representations of GSEA results as ggplot objects for each collection#####

#Store ggplot object (0RvsControl)
#H
HeLa0RvsHeLaControl.GSEAHresult.Dotplot<- ggplot(HeLa0RvsHeLaControl.GSEAHresult.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 30),
        legend.key.size = unit(3, 'cm'),
        legend.text=element_text(size=30))+
  ggtitle("HELA0RvsHELAcontrol GSEA_H Dotplot")

#C2
HeLa0RvsHeLaControl.GSEAC2result.Dotplot<- ggplot(HeLa0RvsHeLaControl.GSEAC2result.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 30),
        legend.key.size = unit(3, 'cm'),
        legend.text=element_text(size=30))+
  ggtitle("HELA0RvsHELAcontrol GSEA_C2 Dotplot")


  #C4
HeLa0RvsHeLaControl.GSEAC4result.Dotplot<- ggplot(HeLa0RvsHeLaControl.GSEAC4result.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 30),
        legend.key.size = unit(3, 'cm'),
        legend.text=element_text(size=30))+
  ggtitle("HELA0RvsHELAcontrol GSEA_C4 Dotplot")

#C5
HeLa0RvsHeLaControl.GSEAC5result.Dotplot<- ggplot(HeLa0RvsHeLaControl.GSEAC5result.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 12),
        legend.key.size = unit(1, 'cm'),
        legend.text=element_text(size=20))+
  ggtitle("HELA0RvsHELAcontrol GSEA_C5 Dotplot")

#C6  
HeLa0RvsHeLaControl.GSEAC6result.Dotplot<- ggplot(HeLa0RvsHeLaControl.GSEAC6result.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 30),
        legend.key.size = unit(3, 'cm'),
        legend.text=element_text(size=30))+
  ggtitle("HELA0RvsHELAcontrol GSEA_C6 Dotplot")



#Store ggplot object (8RvsControl)
HeLa8RvsHeLaControl.GSEAHresult.Dotplot<- ggplot(HeLa8RvsHeLaControl.GSEAHresult.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 30),
        legend.key.size = unit(3, 'cm'),
        legend.text=element_text(size=30))+
  ggtitle("HELA8RvsHELAcontrol GSEA_H Dotplot")

#C2
HeLa8RvsHeLaControl.GSEAC2result.Dotplot<- ggplot(HeLa8RvsHeLaControl.GSEAC2result.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 30),
        legend.key.size = unit(3, 'cm'),
        legend.text=element_text(size=30))+
  ggtitle("HELA8RvsHELAcontrol GSEA_C2 Dotplot")


  #C4
HeLa8RvsHeLaControl.GSEAC4result.Dotplot<- ggplot(HeLa8RvsHeLaControl.GSEAC4result.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 30),
        legend.key.size = unit(3, 'cm'),
        legend.text=element_text(size=30))+
  ggtitle("HELA8RvsHELAcontrol GSEA_C4 Dotplot")

#C5
HeLa8RvsHeLaControl.GSEAC5result.Dotplot<- ggplot(HeLa8RvsHeLaControl.GSEAC5result.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 30),
        legend.key.size = unit(3, 'cm'),
        legend.text=element_text(size=30))+
  ggtitle("HELA8RvsHELAcontrol GSEA_C5 Dotplot")

#C6  
HeLa8RvsHeLaControl.GSEAC6result.Dotplot<- ggplot(HeLa8RvsHeLaControl.GSEAC6result.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 30),
        legend.key.size = unit(3, 'cm'),
        legend.text=element_text(size=30))+
  ggtitle("HELA8RvsHELAcontrol GSEA_C6 Dotplot")

#Store ggplot object (8Rvs0R)
#H
HeLa8RvsHeLa0R.GSEAHresult.Dotplot<- ggplot(HeLa8RvsHeLa0R.GSEAHresult.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 30),
        legend.key.size = unit(3, 'cm'),
        legend.text=element_text(size=30))+
  ggtitle("HeLa8RvsHeLa0R GSEA_H Dotplot")

#C2
HeLa8RvsHeLa0R.GSEAC2result.Dotplot<- ggplot(HeLa8RvsHeLa0R.GSEAC2result.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 30),
        legend.key.size = unit(3, 'cm'),
        legend.text=element_text(size=30))+
  ggtitle("HeLa8RvsHeLa0R GSEA_C2 Dotplot")

  #C4
HeLa8RvsHeLa0R.GSEAC4result.Dotplot<- ggplot(HeLa8RvsHeLa0R.GSEAC4result.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 30),
        legend.key.size = unit(3, 'cm'),
        legend.text=element_text(size=30))+
  ggtitle("HeLa8RvsHeLa0R GSEA_C4 Dotplot")

#C5
HeLa8RvsHeLa0R.GSEAC5result.Dotplot<- ggplot(HeLa8RvsHeLa0R.GSEAC5result.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 30),
        legend.key.size = unit(3, 'cm'),
        legend.text=element_text(size=30))+
  ggtitle("HeLa8RvsHeLa0R GSEA_C5 Dotplot")

#C6  
HeLa8RvsHeLa0R.GSEAC6result.Dotplot<- ggplot(HeLa8RvsHeLa0R.GSEAC6result.tib[1:20,], 
  aes(x=Phenotype, y=ID)) + 
  geom_point(aes(size=setSize, 
                 color =NES , 
                 alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()+ 
  theme(aspect.ratio = 5/1 ,
        text = element_text(size = 30),
        legend.key.size = unit(3, 'cm'),
        legend.text=element_text(size=30))+
  ggtitle("HeLa8RvsHeLa0R GSEA_C6 Dotplot")

#Extract up and down regulated genes defined by LFC cutoffs from the DESeqResults object into a df#####
UP_HELA0RvsHELAcontrol<-filter(as.data.frame(HeLa0RvsHeLaControl), log2FoldChange >=1.5)
DOWN_HELA0RvsHELAcontrol<-filter(as.data.frame(HeLa0RvsHeLaControl), log2FoldChange <=-1.5)
UP_HELA8RvsHELAcontrol<-filter(as.data.frame(HeLa8RvsHeLaControl), log2FoldChange >=1.5)
DOWN_HELA8RvsHELAcontrol<-filter(as.data.frame(HeLa8RvsHeLaControl), log2FoldChange <=-1.5)
UP_HELA8RvsHELA0R<-filter(as.data.frame(HeLa8RvsHeLa0R), log2FoldChange >=1.5)
DOWN_HELA8RvsHELA0R<-filter(as.data.frame(HeLa8RvsHeLa0R), log2FoldChange <=-1.5)

#gprofiler2 enables creation of Manhattan Plots below#####
#Given a list of up/down genes, a model organism and statistical post-hoc correction, returns gprofiler2's GSEA result
HELA0RvsHELAcontrolgost.res_up <- gost(rownames(UP_HELA0RvsHELAcontrol), organism = "hsapiens", correction_method = "fdr")
HELA0RvsHELAcontrolgost.res_down <- gost(rownames(DOWN_HELA0RvsHELAcontrol), organism = "hsapiens", correction_method = "fdr")
HELA8RvsHELAcontrolgost.res_up <- gost(rownames(UP_HELA8RvsHELAcontrol), organism = "hsapiens", correction_method = "fdr")
HELA8RvsHELAcontrolgost.res_down <- gost(rownames(DOWN_HELA8RvsHELAcontrol), organism = "hsapiens", correction_method = "fdr")
HELA8RvsHELA0Rgost.res_up <- gost(rownames(UP_HELA8RvsHELA0R), organism = "hsapiens", correction_method = "fdr")
HELA8RvsHELA0Rgost.res_down <- gost(rownames(DOWN_HELA8RvsHELA0R), organism = "hsapiens", correction_method = "fdr")

#(0RvsControl) Template to use a datatable to filter the NES scores of interest in R. Read in a tibble, NOT the GSEAresult@result 
HeLa0RvsHeLaControl_H_DT<-datatable(HeLa0RvsHeLaControl.GSEAHresult.tib, 
                                     extensions = c('KeyTable', "FixedHeader"), 
                                     caption = 'HeLa0RvsHeLaControl.GSEAHresult',
                                     options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:10), digits=2)

HeLa0RvsHeLaControl_C2_DT<-datatable(HeLa0RvsHeLaControl.GSEAC2result.tib, 
    extensions = c('KeyTable', "FixedHeader"), 
    caption = 'HeLa0RvsHeLaControl.GSEAC2result',
    options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
    formatRound(columns=c(3:10), digits=2)

HeLa0RvsHeLaControl_C4_DT<-datatable(HeLa0RvsHeLaControl.GSEAC4result.tib, 
    extensions = c('KeyTable', "FixedHeader"), 
    caption = 'HeLa0RvsHeLaControl.GSEAC4result',
    options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
    formatRound(columns=c(3:10), digits=2)

HeLa0RvsHeLaControl_C5_DT<-datatable(HeLa0RvsHeLaControl.GSEAC5result.tib, 
    extensions = c('KeyTable', "FixedHeader"), 
    caption = 'HeLa0RvsHeLaControl.GSEAC5result',
    options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
    formatRound(columns=c(3:10), digits=2)

HeLa0RvsHeLaControl_C6_DT<-datatable(HeLa0RvsHeLaControl.GSEAC6result.tib, 
                                     extensions = c('KeyTable', "FixedHeader"), 
                                     caption = 'HeLa0RvsHeLaControl.GSEAC6result',
                                     options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:10), digits=2)

#####Enrichment plots #####

#Enrichment plots(0RvsControl)
#save enrichment plots that mimic the Broad's GSEA results. Export as vector based eps image.
#H NO DATA

#C2
HELA0RvsHELAcontrol.gseaC2.res.eplot<-gseaplot2(
      HeLa0RvsHeLaControl.gseaC2.res, 
      geneSetID = c(1,5,12),
      rel_heights = c(2, 0.5, 0.5),
      base_size = 9,
      pvalue_table = FALSE,
      title = "Top 3 NES Hits: HeLa0RvsHeLaControl GSEA Result C2")

ggsave(file="HeLa0RvsHeLaControl.gseaC2.res.eplotfinaldraft7x4.eps",
       plot = HELA0RvsHELAcontrol.gseaC2.res.eplot,
       device = "eps", 
       units = "in", 
       width= 7,
       height = 4,
       dpi = 300)

#C4
HELA0RvsHELAcontrol.gseaC4.res.eplot<-gseaplot2(
       HeLa0RvsHeLaControl.gseaC4.res, 
       geneSetID = c(3,7,4),
       rel_heights = c(2, 0.5, 0.5),
       base_size = 9,
       pvalue_table = FALSE,
       title = "Top 3 NES Hits: HeLa0RvsHeLaControl GSEA Result C4")

ggsave(file="HeLa0RvsHeLaControl.gseaC4.res.eplotfinaldraft7x4.eps",
       plot = HELA0RvsHELAcontrol.gseaC4.res.eplot,
       device = "eps", 
       units = "in", 
       width= 7,
       height = 4,
       dpi = 300)

#C5
HELA0RvsHELAcontrol.gseaC5.res.eplot<-gseaplot2(
    HeLa0RvsHeLaControl.gseaC5.res, 
    geneSetID = c(3,7,6),
    rel_heights = c(2, 0.5, 0.5),
    base_size = 9,
    pvalue_table = FALSE,
    title = "Top 3 NES Hits: HeLa0RvsHeLaControl GSEA Result C5")

ggsave(file="HeLa0RvsHeLaControl.gseaC5.res.eplotfinaldraft7x4.eps",
       plot = HELA0RvsHELAcontrol.gseaC5.res.eplot,
       device = "eps", 
       units = "in", 
       width= 7,
       height = 4,
       dpi = 300)

#C6
HELA0RvsHELAcontrol.gseaC6.res.eplot<-gseaplot2(
      HeLa0RvsHeLaControl.gseaC6.res, 
      geneSetID = c(1,3,2),
      rel_heights = c(2, 0.5, 0.5),
      base_size = 9,
      pvalue_table = FALSE,
      title = "Top 3 NES Hits: HeLa0RvsHeLaControl GSEA Result C6")

ggsave(file="HELA0RvsHELAcontrol.gseaC6.res.eplotfinaldraft7x4.eps",
       plot = HELA0RvsHELAcontrol.gseaC6.res.eplot,
       device = "eps", 
       units = "in", 
       width= 7,
       height = 4,
       dpi = 300)

#Enrichment plots (8RvsControl)
#save enrichment plots that mimic the Broad's GSEA results. Export as vector based eps image.

#H
HELA8RvsHELAcontrol.gseaH.res.eplot<-gseaplot2(
  HeLa8RvsHeLaControl.gseaH.res, 
  geneSetID = c(1,3,2),
  rel_heights = c(2, 0.5, 0.5),
  base_size = 9,
  pvalue_table = FALSE,
  title = "Top 3 NES Hits: HeLa8RvsHeLaControl GSEA Result H")

ggsave(file="HeLa8RvsHeLaControl.gseaH.res.eplotfinaldraft7x4.eps",
       plot = HELA8RvsHELAcontrol.gseaH.res.eplot,
       device = "eps", 
       units = "in", 
       width= 7,
       height = 4,
       dpi = 300)

#C2
HELA8RvsHELAcontrol.gseaC2.res.eplot<-gseaplot2(
  HeLa8RvsHeLaControl.gseaC2.res, 
  geneSetID = c(1,7,24),
  rel_heights = c(2, 0.5, 0.5),
  base_size = 9,
  pvalue_table = FALSE,
  title = "Top 3 NES Hits: HeLa8RvsHeLaControl GSEAresult C2")

ggsave(file="HeLa8RvsHeLaControl.gseaC2.res.eplotfinaldraft7x4.eps",
       plot = HELA8RvsHELAcontrol.gseaC2.res.eplot,
       device = "eps", 
       units = "in", 
       width= 7,
       height = 4,
       dpi = 300)

#C4
HELA8RvsHELAcontrol.gseaC4.res.eplot<-gseaplot2(
  HeLa8RvsHeLaControl.gseaC4.res, 
  geneSetID = c(4,9,5),
  rel_heights = c(2, 0.5, 0.5),
  base_size = 9,
  pvalue_table = FALSE,
  title = "Top 3 NES Hits: HeLa8RvsHeLaControl GSEAresult C4")

ggsave(file="HELA8RvsHELAcontrol.gseaC4.res.eplotfinaldraft7x4.eps",
       plot = HELA8RvsHELAcontrol.gseaC4.res.eplot,
       device = "eps", 
       units = "in", 
       width= 7,
       height = 4,
       dpi = 300)

#C5
HELA8RvsHELAcontrol.gseaC5.res.eplot<-gseaplot2(
    HeLa8RvsHeLaControl.gseaC5.res, 
    geneSetID = c(5,7,1),
    rel_heights = c(2, 0.5, 0.5),
    base_size = 9,
    pvalue_table = FALSE,
    title = "Top 3 NES Hits: HeLa8RvsHeLaControl GSEAresult C5")

ggsave(file="HeLa8RvsHeLaControl.gseaC5.res.eplotfinaldraft7x4.eps",
       plot = HELA8RvsHELAcontrol.gseaC5.res.eplot,
       device = "eps", 
       units = "in", 
       width= 7,
       height = 4,
       dpi = 300)

#C6
HELA8RvsHELAcontrol.gseaC6.res.eplot<-gseaplot2(
  HeLa8RvsHeLaControl.gseaC6.res, 
  geneSetID = c(1,3,2),
  rel_heights = c(2, 0.5, 0.5),
  base_size = 9,
  pvalue_table = FALSE,
  title = "Top 3 NES Hits: HeLa8RvsHeLaControl GSEAresult C6")

ggsave(file="HeLa8RvsHeLaControl.gseaC6.res.eplotfinaldraft7x4.eps",
       plot = HELA8RvsHELAcontrol.gseaC6.res.eplot,
       device = "eps", 
       units = "in", 
       width= 7,
       height = 4,
       dpi = 300)


#Enrichment plots (8Rvs0R)
#save enrichment plots that mimic the Broad's GSEA results. Export as vector based eps image.

#H
HELA8RvsHELA0R.gseaH.res.eplot<-gseaplot2(
  HeLa8RvsHeLa0R.gseaH.res, 
  geneSetID = c(2,5,1),
  rel_heights = c(2, 0.5, 0.5),
  base_size = 9,
  pvalue_table = FALSE,
  title = "Top 3 NES Hits: HeLa8RvsHeLa0R GSEAresult H")

ggsave(file="HeLa8RvsHeLa0R.gseaH.res.eplotfinaldraft7x4.eps",
       plot = HELA8RvsHELA0R.gseaH.res.eplot,
       device = "eps", 
       units = "in", 
       width= 7,
       height = 4,
       dpi = 300)

#C2
HELA8RvsHELA0R.gseaC2.res.eplot<-gseaplot2(
  HeLa8RvsHeLa0R.gseaC2.res, 
  geneSetID = c(11,29,28),
  rel_heights = c(2, 0.5, 0.5),
  base_size = 9,
  pvalue_table = FALSE,
  title = "Top 3 NES Hits: HeLa8RvsHeLa0R GSEAresult C2")

ggsave(file="HeLa8RvsHeLa0R.gseaC2.res.eplotfinaldraft7x4.eps",
       plot = HELA8RvsHELA0R.gseaC2.res.eplot,
       device = "eps", 
       units = "in", 
       width= 7,
       height = 4,
       dpi = 300)

#C4 NO DATA

#C5
HELA8RvsHELA0R.gseaC5.res.eplot<-gseaplot2(
      HeLa8RvsHeLa0R.gseaC5.res, 
      geneSetID = c(1,16,2),
      rel_heights = c(2, 0.5, 0.5),
      base_size = 9,
      pvalue_table = FALSE,
      title = "Top 3 NES Hits: HeLa8RvsHeLa0R GSEAresult C5")

ggsave(file="HeLa8RvsHeLa0R.gseaC5.res.eplotfinaldraft7x4.eps",
       plot = HELA8RvsHELA0R.gseaC5.res.eplot,
       device = "eps", 
       units = "in", 
       width= 7,
       height = 4,
       dpi = 300)

#C6
HELA8RvsHELA0R.gseaC6.res.eplot<-gseaplot2(
    HeLa8RvsHeLa0R.gseaC6.res, 
    geneSetID = c(3,1,2),
    rel_heights = c(2, 0.5, 0.5),
    base_size = 9,
    pvalue_table = FALSE,
    title = "Top 3 NES Hits: HeLa8RvsHeLa0R GSEAresult C6")

ggsave(file="HeLa8RvsHeLa0R.gseaC6.resfinaldraft7x4.eps",
       plot = HELA8RvsHELA0R.gseaC6.res.eplot,
       device = "eps", 
       units = "in", 
       width= 7,
       height = 4,
       dpi = 300)

#Manhattan plot####
HELA0RvsHELAcontrol_Manhattan_up<-gostplot(HELA0RvsHELAcontrolgost.res_up, interactive = F, capped = T) + labs(title="Manhattan Plot of Up Genes", subtitle = "HeLa0RvsHeLacontrol", caption = "Gene Ontology, KEGG, Reactome, and more")
HELA8RvsHELAcontrol_Manhattan_up<-gostplot(HELA8RvsHELAcontrolgost.res_up, interactive = F, capped = T) + labs(title="Manhattan Plot of Up Genes", subtitle = "HeLa8RvsHeLacontrol", caption = "Gene Ontology, KEGG, Reactome, and more")
HELA8RvsHELA0R_Manhattan_up<-gostplot(HELA8RvsHELA0Rgost.res_up, interactive = F, capped = T) + labs(title="Manhattan Plot of Up Genes", subtitle = "HeLa8RvsHeLa0R", caption = "Gene Ontology, KEGG, Reactome, and more")
HELA0RvsHELAcontrol_Manhattan_down<-gostplot(HELA0RvsHELAcontrolgost.res_down, interactive = F, capped = T) + labs(title="Manhattan Plot of Down Genes", subtitle = "HeLa0RvsHeLacontrol", caption = "Gene Ontology, KEGG, Reactome, and more")
HELA8RvsHELAcontrol_Manhattan_down<-gostplot(HELA8RvsHELAcontrolgost.res_down, interactive = F, capped = T) + labs(title="Manhattan Plot of Down Genes", subtitle = "HeLa8RvsHeLacontrol", caption = "Gene Ontology, KEGG, Reactome, and more")
HELA8RvsHELA0R_Manhattan_down<-gostplot(HELA8RvsHELA0Rgost.res_down, interactive = F, capped = T) + labs(title="Manhattan Plot of Down Genes", subtitle = "HeLa8RvsHeLa0R", caption = "Gene Ontology, KEGG, Reactome, and more")
