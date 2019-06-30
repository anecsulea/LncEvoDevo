##########################################################################

source("parameters.R")

pathDatasets=paste(path, "supplementary_datasets/", sep="")

##########################################################################

splist=c("Mouse", "Rat", "Chicken")

##########################################################################

for(ref in splist){
  print(ref)
  
  for(tg in setdiff(splist, ref)){

    projstats=read.table(paste(pathDatasets, "SupplementaryDataset5/StatsProjections_From",ref,"_To",tg,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

    rownames(projstats)=projstats$GeneID

    save(projstats, file=paste("RData/data.projection.stats.",ref,".",tg,".RData", sep=""))
  }
  
}

##########################################################################




