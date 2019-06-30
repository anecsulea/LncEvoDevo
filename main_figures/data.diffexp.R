#####################################################################################

source("parameters.R")

pathDatasets=paste(path, "supplementary_datasets/", sep="")

#####################################################################################

for(sp in c("Mouse", "Rat", "Chicken")){

  if(sp!="Chicken"){
    de.global.allreads=read.table(paste(pathDatasets, "SupplementaryDataset4/DifferentialExpression_GlobalAgeEffect_AllReads_", sp, ".txt", sep=""), h=T, stringsAsFactors=F)
    rownames(de.global.allreads)=de.global.allreads$GeneID
    
    de.global.resampled=read.table(paste(pathDatasets, "SupplementaryDataset4/DifferentialExpression_GlobalAgeEffect_ResampledReads_", sp, ".txt", sep=""), h=T, stringsAsFactors=F)
    rownames(de.global.resampled)=de.global.resampled$GeneID
  }
  
  de.consecutive.stages=read.table(paste(pathDatasets, "SupplementaryDataset4/DifferentialExpression_ConsecutiveStages_", sp, ".txt", sep=""), h=T, stringsAsFactors=F)
  rownames(de.consecutive.stages)=de.consecutive.stages$GeneID
  
  if(sp!="Chicken"){
    save(list=c("de.global.allreads", "de.global.resampled", "de.consecutive.stages"), file=paste("RData/data.diffexp.",sp,".RData", sep=""))
  } else{
    save(list=c("de.consecutive.stages"), file=paste("RData/data.diffexp.",sp,".RData", sep=""))
  }
}

#####################################################################################



