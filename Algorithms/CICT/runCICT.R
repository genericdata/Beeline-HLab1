library(ppcor)

args <- commandArgs(trailingOnly = T)
inFile <- args[1]
outFile <-  args[2]
gtFile <-  args[3] #ground truth
edgeType <- args[4] ##'ewMIempirical'  'ewMImm' #'Pearson' #'corP'
supervised.positiveClass <- args[5] # c:causal edges, 'c,rc': causal and reversecausal
supervised.negativeClass<- args[6]  # r:random edges, 'rc': reversecausal  
supervised.gtProp<- args[7] #proportion of GT used for classification
supervised.train<- args[8] #proportion of classification records used for training
fp.fn.cost<- args[9]  #False positive to False Negative cost ratio for threshold identification


# input expression data
inputExpr <- read.table(inFile, sep=",", header = 1, row.names = 1)
geneNames <- rownames(inputExpr)
rownames(inputExpr) <- c(geneNames)

# Run CICT 
if(T){
  write.table(data.frame(test=paste('TEST Completed',Sys.time())), outFile, sep = "\t", quote = FALSE, row.names = FALSE)
}

if(F){
  source('CICT')
  CICT_Results=  CICT(x= t(as.matrix(inputExpr)), method = "spearman")
  
  # Write output to a file
  # https://stackoverflow.com/questions/38664241/ranking-and-counting-matrix-elements-in-r
  # Gene1, Gene2, corVal, pValue
  DF = data.frame(Gene1 = CICT_Results$gene1, Gene2 = CICT_Results$gene2
                  , corVal = CICT_Results$estimate, pValue =  CICT_Results$p.value)
  outDF <- DF[order(DF$corVal, decreasing=TRUE), ]
  write.table(outDF, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
}
