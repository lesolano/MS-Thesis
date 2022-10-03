#####load libraries#####
library(topGO)
library("org.Hs.eg.db")
library(lattice)
####import data#####
Genes_HeLa8RvsHeLa0R<-HeLa8RvsHeLa0R@rownames
rawpval_HeLa8RvsHeLa0R<- HeLa8RvsHeLa0R$pvalue
Genelist_HeLa8RvsHeLa0R<- setNames(rawpval_HeLa8RvsHeLa0R, Genes_HeLa8RvsHeLa0R)

GOdata_HeLa8RvsHeLa0R_MF <- new("topGOdata",
                                  ontology = "MF", 
                                  allGenes = Genelist_HeLa8RvsHeLa0R, 
                                  geneSel=function(p) p < 0.01, 
                                  description ="HeLa8RvsHeLa0R_MF",
                                  nodeSize = 10,
                                  annot=annFUN.org, 
                                  mapping="org.Hs.eg.db", 
                                  ID="Ensembl") 
#View GOdata Object
GOdata_HeLa8RvsHeLa0R_MF

#####Run Fisher, Kolmogorov-Smirnov, and custom TopGO elim stat methods #####
resultFisher_HeLa8RvsHeLa0R_MF <- runTest(GOdata_HeLa8RvsHeLa0R_MF, algorithm = "classic", statistic = "fisher")
resultKS_HeLa8RvsHeLa0R_MF <- runTest(GOdata_HeLa8RvsHeLa0R_MF, algorithm = "classic", statistic = "ks")
resultKS.elim_HeLa8RvsHeLa0R_MF <- runTest(GOdata_HeLa8RvsHeLa0R_MF, algorithm = "elim", statistic = "ks")

allRes_HeLa8RvsHeLa0R_MF <- GenTable(GOdata_HeLa8RvsHeLa0R_MF, classicFisher = resultFisher_HeLa8RvsHeLa0R_MF,
                   classicKS = resultKS_HeLa8RvsHeLa0R_MF, elimKS = resultKS.elim_HeLa8RvsHeLa0R_MF,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

pValue.classic_HeLa8RvsHeLa0R_MF <- score(resultKS_HeLa8RvsHeLa0R_MF)
pValue.elim_HeLa8RvsHeLa0R_MF <- score(resultKS.elim_HeLa8RvsHeLa0R_MF)[names(pValue.classic_HeLa8RvsHeLa0R_MF)]
gstat_HeLa8RvsHeLa0R_MF <- termStat(GOdata_HeLa8RvsHeLa0R_MF, names(pValue.classic_HeLa8RvsHeLa0R_MF))
gSize_HeLa8RvsHeLa0R_MF <- gstat_HeLa8RvsHeLa0R_MF$Annotated / max(gstat_HeLa8RvsHeLa0R_MF$Annotated) * 4

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

gCol_HeLa8RvsHeLa0R_MF <- colMap(gstat_HeLa8RvsHeLa0R_MF$Significant)
#####Visualize elim method is more stringent than classic Fisher
plot(pValue.classic_HeLa8RvsHeLa0R_MF, pValue.elim_HeLa8RvsHeLa0R_MF, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize_HeLa8RvsHeLa0R_MF, col = gCol_HeLa8RvsHeLa0R_MF)

sel.go_HeLa8RvsHeLa0R_MF <- names(pValue.classic_HeLa8RvsHeLa0R_MF)[pValue.elim_HeLa8RvsHeLa0R_MF < pValue.classic_HeLa8RvsHeLa0R_MF]
cbind(termStat(GOdata_HeLa8RvsHeLa0R_MF, sel.go_HeLa8RvsHeLa0R_MF),
      elim = pValue.elim_HeLa8RvsHeLa0R_MF[sel.go_HeLa8RvsHeLa0R_MF],
      classic = pValue.classic_HeLa8RvsHeLa0R_MF[sel.go_HeLa8RvsHeLa0R_MF])

showSigOfNodes(GOdata_HeLa8RvsHeLa0R_MF, score(resultKS.elim_HeLa8RvsHeLa0R_MF), firstSigNodes = 3, useInfo = 'all')

printGraph(GOdata_HeLa8RvsHeLa0R_MF, resultKS.elim_HeLa8RvsHeLa0R_MF, firstSigNodes = 3, fn.prefix = "KS.elim_GOdata_HeLa8RvsHeLa0R_MF", useInfo = "all", pdfSW = TRUE)
printGraph(GOdata_HeLa8RvsHeLa0R_MF, resultKS.elim_HeLa8RvsHeLa0R_MF, firstSigNodes = 3, fn.prefix = "KS.elim_GOdata_HeLa8RvsHeLa0R_MF", useInfo = "all", pdfSW = FALSE)
#Export GO summary to CSV
write.csv(allRes_HeLa8RvsHeLa0R_MF, file = "allRes_HeLa8RvsHeLa0R_MF.csv")