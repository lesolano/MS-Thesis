#Script 9

#This script reproduces the differential gene expression process for all HeLa samples with DESeq2

#load libraries####
library(DESeq2)
library(tidyverse)
library('biomaRt')


#import ENSG IDs and read counts Make a table with relevant annotations####
Exp1_countsimport<-read.csv(file = "rawcounts2021.csv", header = TRUE)

#Prepare cts and coldata objects
RawReadCounts<-Exp1_countsimport
rownames(RawReadCounts)<-c(RawReadCounts$gene_id)
cts<-RawReadCounts[,-1]
coldata<-read.csv(file = 'coldata2021.csv',header = TRUE)
row.names(coldata)<-c(coldata$X)
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Condition_NBE)

featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

#factor condition and set levels
dds$Condition_NBE <- factor(dds$Condition_NBE, levels = c("HeLa_Control","HeLa_0R","HeLa_8R",
                                                          "HEK_Control",
                                                          "HEK_0R","HEK_8R",
                                                          "Hep_Control","Hep_0R",
                                                          "Hep_8R"))

#Uncomment to apply prefilter to dds object of readcounts>9
keep <- rowSums(counts(dds)) >= 250
dds <- dds[keep,]
#Run DESEQ function on dds object, and save to the same object
dds <- DESeq(dds)

HeLa0RvsHeLaControl<-results(dds, contrast=c("Condition_NBE","HeLa_0R","HeLa_Control"), alpha = 0.05)
HeLa8RvsHeLaControl <- results(dds, contrast=c("Condition_NBE","HeLa_8R","HeLa_Control"), alpha = 0.05)
HeLa8RvsHeLa0R <- results(dds, contrast=c("Condition_NBE","HeLa_8R","HeLa_0R"), alpha = 0.05)

summary(HeLa0RvsHeLaControl)
summary(HeLa8RvsHeLaControl)
summary(HeLa8RvsHeLa0R)

#Convert ENSG IDs to HGNC IDs for .csv export####
#Define "biomart" to be used
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

#Capture Gene names from DESeq2 result objects
HeLa0RvsHeLaControl_ENSG<-row.names(HeLa0RvsHeLaControl)
HeLa8RvsHeLaControl_ENSG<-row.names(HeLa8RvsHeLaControl)
HeLa8RvsHeLa0R_ENSG<-row.names(HeLa8RvsHeLa0R)

#Create lookup table/retrieve desired annotations.
HeLa0RvsHeLaControl_ENSGlookup <- getBM(
  mart = mart,
  attributes = c('ensembl_gene_id',
                 'gene_biotype','hgnc_symbol'),
  filter = 'ensembl_gene_id',
  values = HeLa0RvsHeLaControl_ENSG,
  uniqueRows = TRUE)

HeLa8RvsHeLaControl_ENSGlookup <- getBM(
  mart = mart,
  attributes = c('ensembl_gene_id',
                 'gene_biotype','hgnc_symbol'),
  filter = 'ensembl_gene_id',
  values = HeLa8RvsHeLaControl_ENSG,
  uniqueRows = TRUE)

HeLa8RvsHeLa0R_ENSGlookup <- getBM(
  mart = mart,
  attributes = c('ensembl_gene_id',
                 'gene_biotype','hgnc_symbol'),
  filter = 'ensembl_gene_id',
  values = HeLa8RvsHeLa0R_ENSG,
  uniqueRows = TRUE)

#prep a column with ENSG IDs from DESeq2 results object and full join merge with lookup table by "ensembl_gene_id"
HeLa0RvsHeLaControl_DEG_colprep<-rownames_to_column(as.data.frame(HeLa0RvsHeLaControl), "ensembl_gene_id")
HeLa0RvsHeLaControl_DEG_Results<-full_join(HeLa0RvsHeLaControl_DEG_colprep,HeLa0RvsHeLaControl_ENSGlookup,by ="ensembl_gene_id", keep=FALSE)

HeLa8RvsHeLaControl_DEG_colprep<-rownames_to_column(as.data.frame(HeLa8RvsHeLaControl), "ensembl_gene_id")
HeLa8RvsHeLaControl_DEG_Results<-full_join(HeLa8RvsHeLaControl_DEG_colprep,HeLa8RvsHeLaControl_ENSGlookup,by ="ensembl_gene_id", keep=FALSE)

HeLa8RvsHeLa0R_DEG_colprep<-rownames_to_column(as.data.frame(HeLa8RvsHeLa0R), "ensembl_gene_id")
HeLa8RvsHeLa0R_DEG_Results<-full_join(HeLa8RvsHeLa0R_DEG_colprep,HeLa8RvsHeLa0R_ENSGlookup,by ="ensembl_gene_id", keep=FALSE)


#Uncomment lines below to export as .CSV
#write.csv(as.data.frame(HeLa0RvsHeLaControl_DEG_Results), file="HeLa0RvsHeLaControl_DEG_Results.csv")
#write.csv(as.data.frame(HeLa8RvsHeLaControl_DEG_Results), file="HeLa8RvsHeLaControl_DEG_Results.csv")
#write.csv(as.data.frame(HeLa8RvsHeLa0R_DEG_Results), file="HeLa8RvsHeLa0R_DEG_Results.csv")

#examples for how to order results by pvalue 
#Create p-value ordered results
HELA0RvsHELAcontrol_Ordered <-HeLa0RvsHeLaControl[order(HeLa0RvsHeLaControl$pvalue),]
HELA8RvsHELAcontrol_Ordered <-HeLa8RvsHeLaControl[order(HeLa8RvsHeLaControl$pvalue),]
HeLa8RvsHeLa0R_Ordered <-HeLa8RvsHeLa0R[order(HeLa8RvsHeLa0R$pvalue),]
