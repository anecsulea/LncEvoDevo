###########################################################################################

source("parameters.R")

pathDatasets=paste(path, "supplementary_datasets/", sep="")

###########################################################################################

alnstat.mr=read.table(paste(pathDatasets, "SupplementaryDataset5/AlignmentStatistics_PredictedOrthoFamilies_Mouse_Rat.txt", sep=""), h=T, stringsAsFactors=F)
rownames(alnstat.mr)=paste(alnstat.mr$ID.Mouse, alnstat.mr$ID.Rat, sep="_")

alnstat.mc=read.table(paste(pathDatasets, "SupplementaryDataset5/AlignmentStatistics_PredictedOrthoFamilies_Mouse_Chicken.txt", sep=""), h=T, stringsAsFactors=F)
rownames(alnstat.mc)=paste(alnstat.mc$ID.Mouse, alnstat.mc$ID.Chicken, sep="_")

alnstat.rc=read.table(paste(pathDatasets, "SupplementaryDataset5/AlignmentStatistics_PredictedOrthoFamilies_Rat_Chicken.txt", sep=""), h=T, stringsAsFactors=F)
rownames(alnstat.rc)=paste(alnstat.rc$ID.Rat, alnstat.rc$ID.Chicken, sep="_")

###########################################################################################

save(list=c("alnstat.mr","alnstat.mc", "alnstat.rc"), file="RData/data.aln.stats.RData")

###########################################################################################
