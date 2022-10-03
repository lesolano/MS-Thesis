#####load libraries#####
library(topGO)
library("org.Hs.eg.db")
library(lattice)
####import data#####
Genes_HeLa8RvsHeLaControl<-HeLa8RvsHeLaControl@rownames
rawpval_HeLa8RvsHeLaControl<- HeLa8RvsHeLaControl$pvalue
Genelist_HeLa8RvsHeLaControl<- setNames(rawpval_HeLa8RvsHeLaControl, Genes_HeLa8RvsHeLaControl)

GOdata_HeLa8RvsHeLaControl_BP <- new("topGOdata",
                                  ontology = "BP", 
                                  allGenes = Genelist_HeLa8RvsHeLaControl, 
                                  geneSel=function(p) p < 0.01, 
                                  description ="HeLa_R0_Genes",
                                  nodeSize = 10,
                                  annot=annFUN.org, 
                                  mapping="org.Hs.eg.db", 
                                  ID="Ensembl") 
#View GOdata Object
GOdata_HeLa8RvsHeLaControl_BP

#####Run Fisher, Kolmogorov-Smirnov, and custom TopGO elim stat methods #####
resultFisher_HeLa8RvsHeLaControl_BP <- runTest(GOdata_HeLa8RvsHeLaControl_BP, algorithm = "classic", statistic = "fisher")
resultKS_HeLa8RvsHeLaControl_BP <- runTest(GOdata_HeLa8RvsHeLaControl_BP, algorithm = "classic", statistic = "ks")
resultKS.elim_HeLa8RvsHeLaControl_BP <- runTest(GOdata_HeLa8RvsHeLaControl_BP, algorithm = "elim", statistic = "ks")

allRes_HeLa8RvsHeLaControl_BP <- GenTable(GOdata_HeLa8RvsHeLaControl_BP, classicFisher = resultFisher_HeLa8RvsHeLaControl_BP,
                   classicKS = resultKS_HeLa8RvsHeLaControl_BP, elimKS = resultKS.elim_HeLa8RvsHeLaControl_BP,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

pValue.classic_HeLa8RvsHeLaControl_BP <- score(resultKS_HeLa8RvsHeLaControl_BP)
pValue.elim_HeLa8RvsHeLaControl_BP <- score(resultKS.elim_HeLa8RvsHeLaControl_BP)[names(pValue.classic_HeLa8RvsHeLaControl_BP)]
gstat_HeLa8RvsHeLaControl_BP <- termStat(GOdata_HeLa8RvsHeLaControl_BP, names(pValue.classic_HeLa8RvsHeLaControl_BP))
gSize_HeLa8RvsHeLaControl_BP <- gstat_HeLa8RvsHeLaControl_BP$Annotated / max(gstat_HeLa8RvsHeLaControl_BP$Annotated) * 4

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

gCol_HeLa8RvsHeLaControl_BP <- colMap(gstat_HeLa8RvsHeLaControl_BP$Significant)
#####Visualize elim method is more stringent than classic Fisher
plot(pValue.classic_HeLa8RvsHeLaControl_BP, pValue.elim_HeLa8RvsHeLaControl_BP, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize_HeLa8RvsHeLaControl_BP, col = gCol_HeLa8RvsHeLaControl_BP)

sel.go_HeLa8RvsHeLaControl_BP <- names(pValue.classic_HeLa8RvsHeLaControl_BP)[pValue.elim_HeLa8RvsHeLaControl_BP < pValue.classic_HeLa8RvsHeLaControl_BP]
cbind(termStat(GOdata_HeLa8RvsHeLaControl_BP, sel.go_HeLa8RvsHeLaControl_BP),
      elim = pValue.elim_HeLa8RvsHeLaControl_BP[sel.go_HeLa8RvsHeLaControl_BP],
      classic = pValue.classic_HeLa8RvsHeLaControl_BP[sel.go_HeLa8RvsHeLaControl_BP])

showSigOfNodes(GOdata_HeLa8RvsHeLaControl_BP, score(resultKS.elim_HeLa8RvsHeLaControl_BP), firstSigNodes = 3, useInfo = 'all')

#uncomment lines below to export pdf or vector based DAG images. Tabular statistical summary available as .csv export


printGraph(GOdata_HeLa8RvsHeLaControl_BP, resultKS.elim_HeLa8RvsHeLaControl_BP, firstSigNodes = 3, fn.prefix = "KS.elim_GOdata_HeLa8RvsHeLaControl_BP", useInfo = "all", pdfSW = TRUE)
printGraph(GOdata_HeLa8RvsHeLaControl_BP, resultKS.elim_HeLa8RvsHeLaControl_BP, firstSigNodes = 3, fn.prefix = "KS.elim_GOdata_HeLa8RvsHeLaControl_BP", useInfo = "all", pdfSW = FALSE)
#Export GO summary to CSV
write.csv(allRes_HeLa8RvsHeLaControl_BP, file = "allRes_HeLa8RvsHeLaControl_BP.csv")