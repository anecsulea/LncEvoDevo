########################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "/supplementary_datasets/", sep="")
pathResults=paste(path, "results/differential_expression/", sep="")

#######################################################################

## R3.5.0
## DESeq2_1.22.0

set.seed(19)

library(DESeq2)

options(stringsAsFactors=F)

########################################################################

ages.order=c("EarlyEmbryo", "LateEmbryo", "Newborn", "Adult", "Aged")

########################################################################

for(sp in c("Mouse", "Rat")){
  print(sp)

 ########################################################################

  geneinfo=read.table(paste(pathDatasets, "SupplementaryDataset1/GeneInfo_",sp,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

  rownames(geneinfo)=geneinfo$GeneID

  mt=geneinfo$GeneID[which(geneinfo$Chr=="MT")] ## remove MT genes before differential expression analysis
  pc=geneinfo$GeneID[which(geneinfo$SelectedProteinCoding=="Yes")]
  lnc=geneinfo$GeneID[which(geneinfo$SelectedLncRNA=="Yes")]

  ## read count
  
  reads=read.table(paste(pathDatasets, "SupplementaryDataset2/UniqueReadCounts_", sp, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

  rownames(reads)=reads$GeneID
  reads=reads[,which(!(colnames(reads)%in%c("GeneID", "ExonicLength")))]

  reads=reads[which(!(rownames(reads)%in%mt)),]

  ########################################################################

  ## tpm values

  tpm=read.table(paste(pathDatasets, "SupplementaryDataset2/KallistoNormalizedTPM_", sp, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(tpm)=tpm$GeneID

  mean.tpm=read.table(paste(pathDatasets, "SupplementaryDataset3/Statistics_AverageTPM_",sp,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(mean.tpm)=mean.tpm$GeneID

  ########################################################################

  samples=colnames(reads)
  tissues=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))
  ages=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[2]))
  ages=unlist(lapply(ages, function(x) substr(x, 1, nchar(x)-1)))

  colnames(reads)=paste("NbReads", colnames(reads), sep=".")
  
  ########################################################################
  
  for(tissue in unique(tissues)){

    print(paste(sp, tissue))
    
    this.reads=reads[,which(tissues==tissue)]
    this.ages=as.factor(ages[which(tissues==tissue)])

    this.sample.order=unlist(lapply(ages.order, function(x) grep(x, samples[which(tissues==tissue)], value=T)))
    this.tissage.order=unique(unlist(lapply(this.sample.order, function(x) substr(x, 1, nchar(x)-1))))
    
    colData=data.frame("condition"=this.ages)
    dds=DESeqDataSetFromMatrix(countData = this.reads, colData=colData, design = ~ condition)
    dds=DESeq(dds,test="LRT", reduced = ~ 1, minReplicatesForReplace=15)

    res=results(dds)
    res=res[order(res$pvalue),]
    res=as.data.frame(res)

    res$GeneID=rownames(res)
    res=cbind(res, this.reads[rownames(res),], mean.tpm[rownames(res), paste("MeanTPM", this.tissage.order, sep=".")])

    colnames(res)[which(colnames(res)=="pvalue")]="PValue"
    colnames(res)[which(colnames(res)=="padj")]="FDR"

    res$GeneType=rep("other", dim(res)[1])
    res$GeneType[which(res$GeneID%in%pc)]="protein_coding"
    res$GeneType[which(res$GeneID%in%lnc)]="lncRNA"
        
    res=res[,c("GeneID", "GeneType", paste("NbReads", this.sample.order, sep="."), grep("MeanTPM", colnames(res), value=T), "PValue", "FDR")]

    write.table(res, file=paste(pathResults, sp, "/DifferentialExpression_AllAges_AllReads_AllGenes_",tissue, ".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  }
}

########################################################################
