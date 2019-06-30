#################################################################################

source("parameters.R")

pathDatasets=paste(path, "supplementary_datasets/", sep="")

##########################################################################

orthofam=read.table(paste(pathDatasets, "SupplementaryDataset5/PredictedOrthoFamilies_Mouse_Rat_Chicken.txt", sep=""), h=T, stringsAsFactors=F) 

##########################################################################

save(orthofam, file="RData/data.ortho.families.RData")

##########################################################################
