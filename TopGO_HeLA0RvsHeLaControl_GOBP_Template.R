#####load libraries#####
library(topGO)
library("org.Hs.eg.db")
library(lattice)
####import data#####
Genes_HeLa0RvsHeLaControl<-HeLa0RvsHeLaControl@rownames
rawpval_HeLa0RvsHeLaControl<- HeLa0RvsHeLaControl$pvalue
Genelist_HeLa0RvsHeLaControl<- setNames(rawpval_HeLa0RvsHeLaControl, Genes_HeLa0RvsHeLaControl)

GOdata_HeLa0RvsHeLaControl_BP <- new("topGOdata",
                                  ontology = "BP", 
                                  allGenes = Genelist_HeLa0RvsHeLaControl, 
                                  geneSel=function(p) p < 0.01, 
                                  description ="HeLa0RvsHeLaControl_BP",
                                  nodeSize = 10,
                                  annot=annFUN.org, 
                                  mapping="org.Hs.eg.db", 
                                  ID="Ensembl") 
#View GOdata Object
GOdata_HeLa0RvsHeLaControl_BP

#####Run Fisher, Kolmogorov-Smirnov, and custom TopGO elim stat methods #####
resultFisher_HeLa0RvsHeLaControl_BP <- runTest(GOdata_HeLa0RvsHeLaControl_BP, algorithm = "classic", statistic = "fisher")
resultKS_HeLa0RvsHeLaControl_BP <- runTest(GOdata_HeLa0RvsHeLaControl_BP, algorithm = "classic", statistic = "ks")
resultKS.elim_HeLa0RvsHeLaControl_BP <- runTest(GOdata_HeLa0RvsHeLaControl_BP, algorithm = "elim", statistic = "ks")

allRes_HeLa0RvsHeLaControl_BP <- GenTable(GOdata_HeLa0RvsHeLaControl_BP, classicFisher = resultFisher_HeLa0RvsHeLaControl_BP,
                   classicKS = resultKS_HeLa0RvsHeLaControl_BP, elimKS = resultKS.elim_HeLa0RvsHeLaControl_BP,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

pValue.classic_HeLa0RvsHeLaControl_BP <- score(resultKS_HeLa0RvsHeLaControl_BP)
pValue.elim_HeLa0RvsHeLaControl_BP <- score(resultKS.elim_HeLa0RvsHeLaControl_BP)[names(pValue.classic_HeLa0RvsHeLaControl_BP)]
gstat_HeLa0RvsHeLaControl_BP <- termStat(GOdata_HeLa0RvsHeLaControl_BP, names(pValue.classic_HeLa0RvsHeLaControl_BP))
gSize_HeLa0RvsHeLaControl_BP <- gstat_HeLa0RvsHeLaControl_BP$Annotated / max(gstat_HeLa0RvsHeLaControl_BP$Annotated) * 4

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

gCol_HeLa0RvsHeLaControl_BP <- colMap(gstat_HeLa0RvsHeLaControl_BP$Significant)
#####Visualize elim method is more stringent than classic Fisher
plot(pValue.classic_HeLa0RvsHeLaControl_BP, pValue.elim_HeLa0RvsHeLaControl_BP, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize_HeLa0RvsHeLaControl_BP, col = gCol_HeLa0RvsHeLaControl_BP)

sel.go_HeLa0RvsHeLaControl_BP <- names(pValue.classic_HeLa0RvsHeLaControl_BP)[pValue.elim_HeLa0RvsHeLaControl_BP < pValue.classic_HeLa0RvsHeLaControl_BP]
cbind(termStat(GOdata_HeLa0RvsHeLaControl_BP, sel.go_HeLa0RvsHeLaControl_BP),
      elim = pValue.elim_HeLa0RvsHeLaControl_BP[sel.go_HeLa0RvsHeLaControl_BP],
      classic = pValue.classic_HeLa0RvsHeLaControl_BP[sel.go_HeLa0RvsHeLaControl_BP])

showSigOfNodes(GOdata_HeLa0RvsHeLaControl_BP, score(resultKS.elim_HeLa0RvsHeLaControl_BP), firstSigNodes = 3, useInfo = 'all')


#uncomment lines below to export pdf or vector based DAG images. Tabular statistical summary available as .csv export

#printGraph(GOdata_HeLa0RvsHeLaControl_BP, resultKS.elim_HeLa0RvsHeLaControl_BP, firstSigNodes = 3, fn.prefix = "KS.elim_GOdata_HeLa0RvsHeLaControl_BP", useInfo = "all", pdfSW = TRUE)
printGraph(GOdata_HeLa0RvsHeLaControl_BP, resultKS.elim_HeLa0RvsHeLaControl_BP, firstSigNodes = 3, fn.prefix = "KS.elim_GOdata_HeLa0RvsHeLaControl_BP", useInfo = "all", pdfSW = FALSE)
#Export GO summary to CSV
write.csv(allRes_HeLa0RvsHeLaControl_BP, file = "allRes_HeLa0RvsHeLaControl_BP.csv")
