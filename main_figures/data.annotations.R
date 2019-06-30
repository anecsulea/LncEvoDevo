####################################################################

source("parameters.R")

pathDatasets=paste(path, "supplementary_datasets/", sep="")

splist=c("Mouse", "Rat", "Chicken")

####################################################################

for(sp in splist){
  print(sp)
  
  allinfo=read.table(paste(pathDatasets,"/SupplementaryDataset1/GeneInfo_",sp,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(allinfo)=allinfo$GeneID
  
  pc=allinfo$GeneID[which(allinfo$SelectedProteinCoding=="Yes")]
  lnc=allinfo$GeneID[which(allinfo$CandidateLncRNA=="Yes")]
    
  save(list=c("pc", "lnc", "allinfo"), file=paste("RData/data.annotations.",sp,".RData", sep=""))
}

####################################################################
