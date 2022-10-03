#####load libraries#####
library(topGO)
library("org.Hs.eg.db")
library(lattice)
####import data#####
Genes_HeLa0RvsHeLaControl<-HeLa0RvsHeLaControl@rownames
rawpval_HeLa0RvsHeLaControl<- HeLa0RvsHeLaControl$pvalue
Genelist_HeLa0RvsHeLaControl<- setNames(rawpval_HeLa0RvsHeLaControl, Genes_HeLa0RvsHeLaControl)

GOdata_HeLa0RvsHeLaControl_MF <- new("topGOdata",
                                  ontology = "MF", 
                                  allGenes = Genelist_HeLa0RvsHeLaControl, 
                                  geneSel=function(p) p < 0.01, 
                                  description ="HeLa0RvsHeLaControl_MF",
                                  nodeSize = 10,
                                  annot=annFUN.org, 
                                  mapping="org.Hs.eg.db", 
                                  ID="Ensembl") 
#View GOdata Object
GOdata_HeLa0RvsHeLaControl_MF

#####Run Fisher, Kolmogorov-Smirnov, and custom TopGO elim stat methods #####
resultFisher_HeLa0RvsHeLaControl_MF <- runTest(GOdata_HeLa0RvsHeLaControl_MF, algorithm = "classic", statistic = "fisher")
resultKS_HeLa0RvsHeLaControl_MF <- runTest(GOdata_HeLa0RvsHeLaControl_MF, algorithm = "classic", statistic = "ks")
resultKS.elim_HeLa0RvsHeLaControl_MF <- runTest(GOdata_HeLa0RvsHeLaControl_MF, algorithm = "elim", statistic = "ks")

allRes_HeLa0RvsHeLaControl_MF <- GenTable(GOdata_HeLa0RvsHeLaControl_MF, classicFisher = resultFisher_HeLa0RvsHeLaControl_MF,
                   classicKS = resultKS_HeLa0RvsHeLaControl_MF, elimKS = resultKS.elim_HeLa0RvsHeLaControl_MF,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

pValue.classic_HeLa0RvsHeLaControl_MF <- score(resultKS_HeLa0RvsHeLaControl_MF)
pValue.elim_HeLa0RvsHeLaControl_MF <- score(resultKS.elim_HeLa0RvsHeLaControl_MF)[names(pValue.classic_HeLa0RvsHeLaControl_MF)]
gstat_HeLa0RvsHeLaControl_MF <- termStat(GOdata_HeLa0RvsHeLaControl_MF, names(pValue.classic_HeLa0RvsHeLaControl_MF))
gSize_HeLa0RvsHeLaControl_MF <- gstat_HeLa0RvsHeLaControl_MF$Annotated / max(gstat_HeLa0RvsHeLaControl_MF$Annotated) * 4

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

gCol_HeLa0RvsHeLaControl_MF <- colMap(gstat_HeLa0RvsHeLaControl_MF$Significant)
#####Visualize elim method is more stringent than classic Fisher
plot(pValue.classic_HeLa0RvsHeLaControl_MF, pValue.elim_HeLa0RvsHeLaControl_MF, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize_HeLa0RvsHeLaControl_MF, col = gCol_HeLa0RvsHeLaControl_MF)

sel.go_HeLa0RvsHeLaControl_MF <- names(pValue.classic_HeLa0RvsHeLaControl_MF)[pValue.elim_HeLa0RvsHeLaControl_MF < pValue.classic_HeLa0RvsHeLaControl_MF]
cbind(termStat(GOdata_HeLa0RvsHeLaControl_MF, sel.go_HeLa0RvsHeLaControl_MF),
      elim = pValue.elim_HeLa0RvsHeLaControl_MF[sel.go_HeLa0RvsHeLaControl_MF],
      classic = pValue.classic_HeLa0RvsHeLaControl_MF[sel.go_HeLa0RvsHeLaControl_MF])

showSigOfNodes(GOdata_HeLa0RvsHeLaControl_MF, score(resultKS.elim_HeLa0RvsHeLaControl_MF), firstSigNodes = 3, useInfo = 'all')

printGraph(GOdata_HeLa0RvsHeLaControl_MF, resultKS.elim_HeLa0RvsHeLaControl_MF, firstSigNodes = 3, fn.prefix = "KS.elim_GOdata_HeLa0RvsHeLaControl_MF", useInfo = "all", pdfSW = TRUE)
printGraph(GOdata_HeLa0RvsHeLaControl_MF, resultKS.elim_HeLa0RvsHeLaControl_MF, firstSigNodes = 3, fn.prefix = "KS.elim_GOdata_HeLa0RvsHeLaControl_MF", useInfo = "all", pdfSW = FALSE)
#Export GO summary to CSV
write.csv(allRes_HeLa0RvsHeLaControl_MF, file = "allRes_HeLa0RvsHeLaControl_MF.csv")