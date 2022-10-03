#####load libraries#####
library(topGO)
library("org.Hs.eg.db")
library(lattice)
####import data#####
Genes_HeLa8RvsHeLaControl<-HeLa8RvsHeLaControl@rownames
rawpval_HeLa8RvsHeLaControl<- HeLa8RvsHeLaControl$pvalue
Genelist_HeLa8RvsHeLaControl<- setNames(rawpval_HeLa8RvsHeLaControl, Genes_HeLa8RvsHeLaControl)

GOdata_HeLa8RvsHeLaControl_CC <- new("topGOdata",
                                  ontology = "CC", 
                                  allGenes = Genelist_HeLa8RvsHeLaControl, 
                                  geneSel=function(p) p < 0.01, 
                                  description ="HeLa8RvsHeLaControl_CC",
                                  nodeSize = 10,
                                  annot=annFUN.org, 
                                  mapping="org.Hs.eg.db", 
                                  ID="Ensembl") 
#View GOdata Object
GOdata_HeLa8RvsHeLaControl_CC

#####Run Fisher, Kolmogorov-Smirnov, and custom TopGO elim stat methods #####
resultFisher_HeLa8RvsHeLaControl_CC <- runTest(GOdata_HeLa8RvsHeLaControl_CC, algorithm = "classic", statistic = "fisher")
resultKS_HeLa8RvsHeLaControl_CC <- runTest(GOdata_HeLa8RvsHeLaControl_CC, algorithm = "classic", statistic = "ks")
resultKS.elim_HeLa8RvsHeLaControl_CC <- runTest(GOdata_HeLa8RvsHeLaControl_CC, algorithm = "elim", statistic = "ks")

allRes_HeLa8RvsHeLaControl_CC <- GenTable(GOdata_HeLa8RvsHeLaControl_CC, classicFisher = resultFisher_HeLa8RvsHeLaControl_CC,
                   classicKS = resultKS_HeLa8RvsHeLaControl_CC, elimKS = resultKS.elim_HeLa8RvsHeLaControl_CC,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

pValue.classic_HeLa8RvsHeLaControl_CC <- score(resultKS_HeLa8RvsHeLaControl_CC)
pValue.elim_HeLa8RvsHeLaControl_CC <- score(resultKS.elim_HeLa8RvsHeLaControl_CC)[names(pValue.classic_HeLa8RvsHeLaControl_CC)]
gstat_HeLa8RvsHeLaControl_CC <- termStat(GOdata_HeLa8RvsHeLaControl_CC, names(pValue.classic_HeLa8RvsHeLaControl_CC))
gSize_HeLa8RvsHeLaControl_CC <- gstat_HeLa8RvsHeLaControl_CC$Annotated / max(gstat_HeLa8RvsHeLaControl_CC$Annotated) * 4

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

gCol_HeLa8RvsHeLaControl_CC <- colMap(gstat_HeLa8RvsHeLaControl_CC$Significant)
#####Visualize elim method is more stringent than classic Fisher
plot(pValue.classic_HeLa8RvsHeLaControl_CC, pValue.elim_HeLa8RvsHeLaControl_CC, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize_HeLa8RvsHeLaControl_CC, col = gCol_HeLa8RvsHeLaControl_CC)

sel.go_HeLa8RvsHeLaControl_CC <- names(pValue.classic_HeLa8RvsHeLaControl_CC)[pValue.elim_HeLa8RvsHeLaControl_CC < pValue.classic_HeLa8RvsHeLaControl_CC]
cbind(termStat(GOdata_HeLa8RvsHeLaControl_CC, sel.go_HeLa8RvsHeLaControl_CC),
      elim = pValue.elim_HeLa8RvsHeLaControl_CC[sel.go_HeLa8RvsHeLaControl_CC],
      classic = pValue.classic_HeLa8RvsHeLaControl_CC[sel.go_HeLa8RvsHeLaControl_CC])

#THIS IS BUGGED!!!! 
showSigOfNodes(GOdata_HeLa8RvsHeLaControl_CC, score(resultKS.elim_HeLa8RvsHeLaControl_CC), firstSigNodes = 3, useInfo = 'all')

printGraph(GOdata_HeLa8RvsHeLaControl_CC, resultKS.elim_HeLa8RvsHeLaControl_CC, firstSigNodes = 3, fn.prefix = "KS.elim_GOdata_HeLa8RvsHeLaControl_CC", useInfo = "all", pdfSW = TRUE)
printGraph(GOdata_HeLa8RvsHeLaControl_CC, resultKS.elim_HeLa8RvsHeLaControl_CC, firstSigNodes = 3, fn.prefix = "KS.elim_GOdata_HeLa8RvsHeLaControl_CC", useInfo = "all", pdfSW = FALSE)
#Export GO summary to CSV
write.csv(allRes_HeLa8RvsHeLaControl_CC, file = "allRes_HeLa8RvsHeLaControl_CC.csv")