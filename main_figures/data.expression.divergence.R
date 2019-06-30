#########################################################################

source("parameters.R")

pathDatasets=paste(path, "supplementary_datasets/", sep="")

#########################################################################

expdiv=read.table(paste(pathDatasets, "SupplementaryDataset7/ExpressionDivergence_OrthoGenes_Mouse_Rat.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

save(expdiv, file="RData/data.expression.divergence.RData")

########################################################################

