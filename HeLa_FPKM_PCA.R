#load Libraries
library(tidyverse)
library(ggfortify)

#load normalized counts from E1 and E2
HEKHELAHEP_E1FPKM<-read.csv("E1v2.csv", header = TRUE)
HEKHELAHEP_E2FPKM<-read.csv("E2v2.csv", header = TRUE)

#Assign row names
rownames(HEKHELAHEP_E1FPKM)<-HEKHELAHEP_E1FPKM$geneID
rownames(HEKHELAHEP_E2FPKM)<-HEKHELAHEP_E2FPKM$geneID

#full join by ENSGID and save gene ENSG IDs
HEKHELAHEP_E1E2<-full_join(HEKHELAHEP_E1FPKM, HEKHELAHEP_E2FPKM, by ='geneID')
geneID<-HEKHELAHEP_E1E2$geneID

#grep out HeLa column names from the full join and join with gene IDs
HELAcols<-HEKHELAHEP_E1E2[,grepl("HeLa",names(HEKHELAHEP_E1E2))]
HELA_E1E2<-cbind(geneID,HELAcols)

#Transpose and create a dataframe. Need Transcript IDs as row names and Samples as column names
HELA_E1E2t<-as.data.frame(t(HELA_E1E2))

#Assign Gene IDs (first column entires of the original read-in table) to dataframe transpose columns 
colnames(HELA_E1E2t) <- HELA_E1E2[,1]

#Remove the first row (containing the ENSG IDs) that was saved earlier to rownames)
HELA_E1E2t<-HELA_E1E2t[c(-1),]

HELAsampleconds<-read.csv2("HELAsampleconditions.csv",header = TRUE, sep = ",")
Experiment<-HELAsampleconds$Experiment
Condition<-HELAsampleconds$Condition
ConditionbyExp<-HELAsampleconds$ConditionbyExp

HELA_E1E2t["Condition"]<- Condition

HELA_E1E2t["Experiment"]<- Experiment

HELA_E1E2t["ConditionbyExp"]<- ConditionbyExp



#HEKHELAHEP_E1E2t_poscols = HEKHELAHEP_E1E2t[,colSums(HEKHELAHEP_E1E2t) > 0]
HELA_E1E2t_posvar<-HELA_E1E2t[,which(apply(HELA_E1E2t, 2, var) != 0)]
#need to convert characters in FPKMt into numbers
HELA_E1E2_pcainputFPKM<-data.matrix(HELA_E1E2t_posvar)

#Run PCA with prcomp. BreastLiverKidney_pca_res should be an S4 object with stats, loadings etc
HELA_pca_res <- prcomp(HELA_E1E2_pcainputFPKM, scale. = TRUE)

#plots the PCA SCORE/SCREE####

#No knowledge of sample type/conditions
blindPCA<-autoplot(HELA_pca_res)

#Color by Condition or batch or both
HELAConditionPCA<-autoplot(HELA_pca_res , data = HELA_E1E2t, colour = 'Condition') + theme(text = element_text(size = 40))
HELAbyExpPCA<-autoplot(HELA_pca_res , data = HELA_E1E2t, colour = 'Experiment')
HELAConditionbyExpPCA<-autoplot(HELA_pca_res , data = HELA_E1E2t, colour = 'ConditionbyExp')

#Scree Plot
#proportion of variance in $importance[2, ] and cumulative variance in $importance[3, ]
HELAsummaryObj<-summary(HELA_pca_res)

#individual contributions
HELA_screeplot<-qplot(c(1:25), HELAsummaryObj$importance[2, ]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 0.2) +
  xlim(1,6) +
  theme(text = element_text(size = 40))  

#cumulative contributions
HELA_cumulative_screeplot<-qplot(c(1:25), HELAsummaryObj$importance[3, ]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Cumulative Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1) +
  xlim(1,6) +
  theme(text = element_text(size = 40))  

#Combined plot
HELA_combined_screeplot<-qplot(c(1:25), HELAsummaryObj, aes(y=HELAsummaryObj$importance[3, ])) + 
  geom_line(aes(y=HELAsummaryObj$importance[2, ])) + 
  geom_line(aes(y=HELAsummaryObj$importance[3, ])) +
  geom_point(aes(y=HELAsummaryObj$importance[2, ])) +
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("HeLa Scree Plot and Cumulative Variance Explained") +
  ylim(0, 0.5) +
  xlim(1,6) +
  theme(text = element_text(size = 40))  

#Biplot works but genes are eigenvectors and it's too many arrows to be useful. 
#Might be useful if only top X genes were plotted, but score plots already gives us PCXvsPCY.

#eigenvectors are negative in R by default. Multiply by -1 to reverse
#rotation contains gene contributions to PCs, X contains sample contributions
#HELA_pca_res$rotation <- -1*HELA_pca_res$rotation
#HELA_pca_res$x <- -1*HELA_pca_res$x

#biplot(HELA_pca_res, scale =0)
