library(AUCell)
library(GSEABase)
library(data.table)

# importing the data file in the format gene id, gene name, gene expression per cell
geoData <- fread("GSE147405_MCF7_TGFB1.txt", sep="\t")
#geoData <- subset(geoData,select = -c(Gene_Name))
geneNames <- toupper(unname(unlist(geoData[,1, with=FALSE])))
exprMatrix <- as.matrix(geoData[,-1, with=FALSE])
dim(exprMatrix)
rownames(exprMatrix) <- geneNames

geneSets <- getGmt("./../signatures/signatures.gmt")
geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix))
cbind(nGenes(geneSets))
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep=""))
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, nCores=1)
x <- getAUC(cells_AUC)[,]
write.table(x,sep = "\t",file = "GSE147405_MCF7_TGFB1.tsv",quote = FALSE)