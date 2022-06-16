args <- commandArgs(trailingOnly = T)
inFile <- args[1]
outFile <-  args[2]

# input expression data
inputExpr <- read.table(inFile, sep=",", header = 1, row.names = 1)
geneNames <- rownames(inputExpr)
rownames(inputExpr) <- c(geneNames)

set.seed(runif(1,1,10000))
DF = expand.grid(Gene1=geneNames,Gene2=geneNames)
DF$predVal = runif(nrow(DF))
outDF = DF[order(DF$predVal, decreasing=TRUE), ]
write.table(outDF, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
