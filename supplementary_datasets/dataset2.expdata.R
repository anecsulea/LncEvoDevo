##########################################################################

path="LncEvoDevo/"
pathExpression=paste(path, "results/expression_estimation/", sep="")
pathResampledReads=paste(path, "results/read_counts_resampled_noMT/", sep="")
pathResults=paste(path, "supplementary_datasets/SupplementaryDataset2/", sep="")

options(scipen=999) ## remove scientific notation ## options(scipen=0) to get it back

nbresampled=c("60M","60M", "50M")
names(nbresampled)=c("Mouse", "Rat", "Chicken")

##########################################################################

for(sp in c("Mouse", "Rat", "Chicken")){
  print(sp)

  ##################################
  
  counts.subread=read.table(paste(pathExpression, sp, "/AllSamples_UniqueReadCounts_StringTie_MainStrain.txt", sep=""), h=T, stringsAsFactors=F)
  counts.subread=counts.subread[, -which(colnames(counts.subread)=="ExonicLength")]
  
  counts.kallisto=read.table(paste(pathExpression, sp, "/AllSamples_KallistoEstimatedCounts_StringTie_MainStrain.txt", sep=""), h=T, stringsAsFactors=F)
  
  ##################################

  resampled=read.table(paste(pathResampledReads, sp, "/AllSamples_UniqueReadCounts_StringTie.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  geneid=resampled$GeneID
  resampled=resampled[,-1]
  tissage=as.factor(unlist(lapply(colnames(resampled), function(x) substr(x,1, nchar(x)-1))))
  resampled.tissage=t(apply(resampled, 1, function(x) tapply(as.numeric(x), tissage, sum)))
  colnames(resampled.tissage)=levels(tissage)

  resampled.tissage=data.frame("GeneID"=geneid, resampled.tissage, stringsAsFactors=F)

  ##################################

  raw.tpm=read.table(paste(pathExpression, sp, "/AllSamples_KallistoTPM_StringTie_MainStrain.txt", sep=""), h=T, stringsAsFactors=F)
    
  norm.tpm=read.table(paste(pathExpression, sp, "/AllSamples_KallistoNormalizedTPM_StringTie_MainStrain.txt", sep=""), h=T, stringsAsFactors=F)
    
  ##################################

  ### output

  write.table(counts.subread, file=paste(pathResults, "UniqueReadCounts_",sp,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

  write.table(counts.kallisto, file=paste(pathResults, "KallistoEstimatedCounts_",sp,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  
  write.table(raw.tpm, file=paste(pathResults, "KallistoRawPM_",sp,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  
  write.table(norm.tpm, file=paste(pathResults, "KallistoNormalizedTPM_",sp,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

  write.table(resampled.tissage, file=paste(pathResults, "UniqueReadCounts_Downsampled",nbresampled[sp],"_",sp,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

  ##################################
}

#########################################################################

  
