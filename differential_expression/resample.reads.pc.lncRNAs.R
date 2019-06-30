########################################################################

path="LncEvoDevo/"
pathResults=paste(path, "results/differential_expression/", sep="")
pathDatasets=paste(path, "supplementary_datasets/", sep="")

########################################################################################

set.seed(19)

ages.order=c("EarlyEmbryo", "LateEmbryo", "Newborn", "Adult", "Aged")

splist=c("Mouse", "Rat", "Chicken")

########################################################################################

for(sp in splist){
  print(sp)

  ## gene info
  
  geneinfo=read.table(paste(pathDatasets, "SupplementaryDataset1/GeneInfo_",sp,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

  rownames(geneinfo)=geneinfo$GeneID

  mt=geneinfo$GeneID[which(geneinfo$Chr=="MT")] ## remove MT genes before differential expression analysis
  pc=geneinfo$GeneID[which(geneinfo$SelectedProteinCoding=="Yes" & geneinfo$Chr!="MT")]
  lnc=geneinfo$GeneID[which(geneinfo$SelectedLncRNA=="Yes" & geneinfo$Chr!="MT")]

  ## read counts
  
  reads=read.table(paste(pathDatasets, "SupplementaryDataset2/UniqueReadCounts_", sp, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

  rownames(reads)=reads$GeneID
  reads=reads[,which(!(colnames(reads)%in%c("GeneID")))]

  reads=reads[which(!(rownames(reads)%in%mt)),] ## remove mitochondrial genes

  ## sample info
  
  samples=colnames(reads)
  tissues=unlist(lapply(samples, function(x) unlist(strsplit(x,split="_"))[1]))
  ages=unlist(lapply(samples, function(x) unlist(strsplit(x,split="_"))[2]))
  ages=unlist(lapply(ages, function(x) substr(x,1,nchar(x)-1)))

  ########################################################################################

  reads.pc=reads[which(rownames(reads)%in%pc),]
  reads.lnc=reads[which(rownames(reads)%in%lnc),]

  nbpc=dim(reads.pc)[1]
  nblnc=dim(reads.lnc)[1]

  reads.pc.resampled=reads.pc

  for(s in samples){
    print(s)
    
    mean.lnc=sum(reads.lnc[,s])/dim(reads.lnc)[1]
    mean.pc=sum(reads.pc[,s])/dim(reads.pc)[1]

    tot.pc=sum(reads.pc[,s])
    prop.pc=reads.pc[,s]/sum(reads.pc[,s])

    newtot=floor(tot.pc*mean.lnc/mean.pc)
  
    newreads=sample(1:nbpc, size=newtot,  prob=prop.pc, replace=TRUE)
    newcounts=as.numeric(table(factor(newreads, levels=1:nbpc)))

    reads.pc.resampled[,s]=newcounts
  }
  
  save(list=c("reads.lnc", "reads.pc", "reads.pc.resampled", "samples", "tissues", "ages"), file=paste("RData/reads.resampling.pc.lncRNA.",sp,".RData",sep=""))
}

########################################################################
