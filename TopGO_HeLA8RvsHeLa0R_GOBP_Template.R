#####load libraries#####
library(topGO)
library("org.Hs.eg.db")
library(lattice)
####import data#####
Genes_HeLa8RvsHeLa0R<-HeLa8RvsHeLa0R@rownames
rawpval_HeLa8RvsHeLa0R<- HeLa8RvsHeLa0R$pvalue
Genelist_HeLa8RvsHeLa0R<- setNames(rawpval_HeLa8RvsHeLa0R, Genes_HeLa8RvsHeLa0R)

GOdata_HeLa8RvsHeLa0R_BP <- new("topGOdata",
                                  ontology = "BP", 
                                  allGenes = Genelist_HeLa8RvsHeLa0R, 
                                  geneSel=function(p) p < 0.01, 
                                  description ="HeLa8RvsHeLa0R_BP",
                                  nodeSize = 10,
                                  annot=annFUN.org, 
                                  mapping="org.Hs.eg.db", 
                                  ID="Ensembl") 
#View GOdata Object
GOdata_HeLa8RvsHeLa0R_BP

#####Run Fisher, Kolmogorov-Smirnov, and custom TopGO elim stat methods #####
resultFisher_HeLa8RvsHeLa0R_BP <- runTest(GOdata_HeLa8RvsHeLa0R_BP, algorithm = "classic", statistic = "fisher")
resultKS_HeLa8RvsHeLa0R_BP <- runTest(GOdata_HeLa8RvsHeLa0R_BP, algorithm = "classic", statistic = "ks")
resultKS.elim_HeLa8RvsHeLa0R_BP <- runTest(GOdata_HeLa8RvsHeLa0R_BP, algorithm = "elim", statistic = "ks")

allRes_HeLa8RvsHeLa0R_BP <- GenTable(GOdata_HeLa8RvsHeLa0R_BP, classicFisher = resultFisher_HeLa8RvsHeLa0R_BP,
                   classicKS = resultKS_HeLa8RvsHeLa0R_BP, elimKS = resultKS.elim_HeLa8RvsHeLa0R_BP,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

pValue.classic_HeLa8RvsHeLa0R_BP <- score(resultKS_HeLa8RvsHeLa0R_BP)
pValue.elim_HeLa8RvsHeLa0R_BP <- score(resultKS.elim_HeLa8RvsHeLa0R_BP)[names(pValue.classic_HeLa8RvsHeLa0R_BP)]
gstat_HeLa8RvsHeLa0R_BP <- termStat(GOdata_HeLa8RvsHeLa0R_BP, names(pValue.classic_HeLa8RvsHeLa0R_BP))
gSize_HeLa8RvsHeLa0R_BP <- gstat_HeLa8RvsHeLa0R_BP$Annotated / max(gstat_HeLa8RvsHeLa0R_BP$Annotated) * 4

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

gCol_HeLa8RvsHeLa0R_BP <- colMap(gstat_HeLa8RvsHeLa0R_BP$Significant)
#####Visualize elim method is more stringent than classic Fisher
plot(pValue.classic_HeLa8RvsHeLa0R_BP, pValue.elim_HeLa8RvsHeLa0R_BP, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize_HeLa8RvsHeLa0R_BP, col = gCol_HeLa8RvsHeLa0R_BP)

sel.go_HeLa8RvsHeLa0R_BP <- names(pValue.classic_HeLa8RvsHeLa0R_BP)[pValue.elim_HeLa8RvsHeLa0R_BP < pValue.classic_HeLa8RvsHeLa0R_BP]
cbind(termStat(GOdata_HeLa8RvsHeLa0R_BP, sel.go_HeLa8RvsHeLa0R_BP),
      elim = pValue.elim_HeLa8RvsHeLa0R_BP[sel.go_HeLa8RvsHeLa0R_BP],
      classic = pValue.classic_HeLa8RvsHeLa0R_BP[sel.go_HeLa8RvsHeLa0R_BP])

showSigOfNodes(GOdata_HeLa8RvsHeLa0R_BP, score(resultKS.elim_HeLa8RvsHeLa0R_BP), firstSigNodes = 3, useInfo = 'all')

#uncomment lines below to export pdf or vector based DAG images. Tabular statistical summary available as .csv export


printGraph(GOdata_HeLa8RvsHeLa0R_BP, resultKS.elim_HeLa8RvsHeLa0R_BP, firstSigNodes = 3, fn.prefix = "KS.elim_GOdata_HeLa8RvsHeLa0R_BP", useInfo = "all", pdfSW = TRUE)
printGraph(GOdata_HeLa8RvsHeLa0R_BP, resultKS.elim_HeLa8RvsHeLa0R_BP, firstSigNodes = 3, fn.prefix = "KS.elim_GOdata_HeLa8RvsHeLa0R_BP", useInfo = "all", pdfSW = FALSE)
#Export GO summary to CSV
write.csv(allRes_HeLa8RvsHeLa0R_BP, file = "allRes_HeLa8RvsHeLa0R_BP.csv")