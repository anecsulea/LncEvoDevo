####################################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/",sep="")
pathTables=paste(path, "supplementary_tables/",sep="")

####################################################################################

species=c()
sample=c()
nbpc=c()
nblnc=c()

minTPM=1

tissue.order=c("Brain","Kidney","Liver","Testis")
age.order=c("EarlyEmbryo", "LateEmbryo", "Newborn", "Adult", "Aged")
newage=c("MidStageEmbryo","LateEmbryo", "Newborn", "Adult", "Aged")

names(newage)=age.order

for(sp in c("Mouse", "Rat", "Chicken")){
  info=read.table(paste(pathDatasets, "SupplementaryDataset1/GeneInfo_",sp,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  pc=info$GeneID[which(info$SelectedProteinCoding=="Yes")]
  lnc=info$GeneID[which(info$CandidateLncRNA=="Yes")]
  
  tpm=read.table(paste(pathDatasets, "SupplementaryDataset3/Statistics_AverageTPM_",sp,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  rownames(tpm)=tpm$GeneID

  for(tiss in tissue.order){
    for(age in age.order){
      this.nbpc=length(which(tpm[pc,paste("MeanTPM.",tiss,"_",age,sep="")]>=minTPM))
      this.nblnc=length(which(tpm[lnc,paste("MeanTPM.",tiss,"_",age,sep="")]>=minTPM))

      nbpc=c(nbpc, this.nbpc)
      nblnc=c(nblnc, this.nblnc)
      species=c(species, sp)
      sample=c(sample, paste(tiss, newage[age], sep="_"))
      
    } 
  }
}

results=data.frame("Species"=species, "Sample"=sample, "NbProteinCoding"=nbpc, "NbLncRNAs"=nblnc, stringsAsFactors=F)
results=results[which(results$NbProteinCoding>0),]

write.table(results, file=paste(pathTables, "SupplementaryTable5.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

####################################################################################


