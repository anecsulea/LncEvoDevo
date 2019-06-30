############################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathResults=paste(path, "results/mutual_information_network/", sep="")

options(scipen=999)

maxFDR=0.001

############################################################################

for(type in c("lncRNAs_only")){#, "all_genes"
  mouse=read.table(paste(pathResults, type, "/Mouse/filtered_network_FDR",maxFDR,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  mouse$ID=apply(mouse[,c("Regulator", "Target")], 1, function(x) paste(sort(as.character(x)), collapse="-"))
  dupli=which(duplicated(mouse$ID))
  
  if(length(dupli)>0){
    mouse=mouse[-dupli,]
  }

  rownames(mouse)=mouse$ID
  
  ## rat
  
  rat=read.table(paste(pathResults, type, "/Rat/filtered_network_FDR",maxFDR,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rat$ID=apply(rat[,c("Regulator", "Target")], 1, function(x) paste(sort(as.character(x)),  collapse="-"))
  dupli=which(duplicated(rat$ID))
  
  if(length(dupli)>0){
    rat=rat[-dupli,]
  }

  rownames(rat)=rat$ID
  
  ## intersect
  common=intersect(mouse$ID, rat$ID)

  results=data.frame("ID"=common, "ID1"=mouse[common,"Regulator"], "ID2"=mouse[common, "Target"], "MI.Mouse"=mouse[common,"MI"], "MI.Rat"=rat[common, "MI"], "PValue.Mouse"=mouse[common, "pvalue"], "PValue.Rat"=rat[common, "pvalue"], "FDR.Mouse"=mouse[common, "FDR"], "FDR.Rat"=rat[common, "FDR"], stringsAsFactors=F) 

  write.table(results, file=paste(pathResults, type, "/common_interactions_FDR",maxFDR,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

}

############################################################################
