#This script reproduces the differential gene expression process for all HeLa samples with DESeq2

#load libraries####
library(DESeq2)
library(tidyverse)

#import ENSG IDs and raw read counts####
Exp1_countsimport<-read.csv(file = "rawcounts2021.csv", header = TRUE)
Exp2_countsimport<-read.csv(file = "rawcounts2022.csv", header = TRUE)

#Create a single dataframe from the raw counts data with shared ENSG ID column
Exp1_2_countsimport<-full_join(Exp1_countsimport,Exp2_countsimport,by ="gene_id", keep=FALSE)

#Define cts and coldata objects with appropriate annotations 
RawReadCounts<-Exp1_2_countsimport
rownames(RawReadCounts)<-c(RawReadCounts$gene_id)
cts<-RawReadCounts[,-1]
coldata<-read.csv(file = 'coldata220924.csv',header = TRUE)
row.names(coldata)<-c(coldata$X)
all(rownames(coldata) %in% colnames(cts))
#If this line returns FALSE, uncomment the indicated lines of code below
all(rownames(coldata) == colnames(cts))

#This ensures coldata rownames match the column names of the cts object exactly
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

#Create dds object and add metadata
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Condition_NBE)

featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

#factor condition and set levels
dds$Condition_NBE <- factor(dds$Condition_NBE, levels = c("HeLa_0R","HeLa_8R",
                                                          "HeLa_Control","HEK_Control",
                                                          "HEK_0R","HEK_8R",
                                                          "Hep_Control","Hep_0R",
                                                          "Hep_8R"))

#apply prefilter to dds object of readcounts>9
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#RunDeseq2 by contrast, store results to descriptively named results objects
#Use contrast to designate comparisons of interest.  
dds <- DESeq(dds)
HeLa0RvsHeLaControl<-results(dds, contrast=c("Condition_NBE","HeLa_0R","HeLa_Control"), alpha = 0.05)
HeLa8RvsHeLaControl <- results(dds, contrast=c("Condition_NBE","HeLa_8R","HeLa_Control"), alpha = 0.05)
HeLa8RvsHeLa0R <- results(dds, contrast=c("Condition_NBE","HeLa_8R","HeLa_0R"), alpha = 0.05)

#Uncomment lines below to export DEG results to a .CSV
#write.csv(as.data.frame(HeLa0RvsHeLaControl), file="HeLa0RvsHeLaControl_NovogeneRawCounts.csv")
#write.csv(as.data.frame(HeLa8RvsHeLaControl), file="HeLa8RvsHeLaControl_NovogeneRawCounts.csv")
#write.csv(as.data.frame(HeLa8RvsHeLa0R), file="HeLa8RvsHeLa0R_NovogeneRawCounts.csv")

#examples for how to order results by pvalue 
#Create p-value ordered results
HELA0RvsHELAcontrol_Ordered <-HeLa0RvsHeLaControl[order(HeLa0RvsHeLaControl$pvalue),]
HELA8RvsHELAcontrol_Ordered <-HeLa8RvsHeLaControl[order(HeLa8RvsHeLaControl$pvalue),]
HeLa8RvsHeLa0R_Ordered <-HeLa8RvsHeLa0R[order(HeLa8RvsHeLa0R$pvalue),]
