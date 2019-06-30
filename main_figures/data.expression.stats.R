####################################################################

source("parameters.R")

pathDatasets=paste(path, "supplementary_datasets/", sep="")

splist=c("Mouse", "Rat", "Chicken")

####################################################################

for(sp in splist){
  stats=read.table(paste(pathDatasets, "SupplementaryDataset3/Statistics_AverageTPM_",sp,".txt", sep=""), h=T, stringsAsFactors=F) 
  rownames(stats)=stats$GeneID
  stats=stats[,-1]

  save(stats, file=paste("RData/data.tpm.stats.",sp,".RData", sep=""))
}

####################################################################
