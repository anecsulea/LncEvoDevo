########################################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathDiffExp=paste(path, "results/differential_expression/", sep="")

#options(scipen=999) ## remove scientific notation ## options(scipen=0) to get it back

age.order=c("EarlyEmbryo", "LateEmbryo", "Newborn", "Adult", "Aged")

########################################################################################

## differential expression only for selected protein-coding genes and lncRNAs

########################################################################################

for(sp in c("Mouse", "Rat", "Chicken")){
  geneinfo=read.table(paste(pathDatasets, "SupplementaryDataset1/GeneInfo_",sp,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(geneinfo)=geneinfo$GeneID

  pc=geneinfo$GeneID[which(geneinfo$SelectedProteinCoding=="Yes")]
  lnc=geneinfo$GeneID[which(geneinfo$SelectedLncRNA=="Yes")]
  allgenes=c(pc, lnc)
  genetype=c(rep("ProteinCoding", length(pc)), rep("LncRNA", length(lnc)))

  de.global.allreads=list()
  de.global.resampled=list()
  de.consecutive.stages=list()
  
  for(tissue in c("Brain", "Kidney", "Liver", "Testis")){
    print(paste(sp, tissue))
    
    ## global difference between stages, all reads
    
    pathde=paste(pathDiffExp, sp, "/DifferentialExpression_AllAges_AllReads_PCLncRNAs_",tissue,".txt", sep="")
    
    if(file.exists(pathde)){
      this.de.global=read.table(pathde, h=T, stringsAsFactors=F, sep="\t", quote="\"")
      rownames(this.de.global)=this.de.global$GeneID
      
      de.global.allreads[[paste("FDR", tissue, sep=".")]]=this.de.global[allgenes, "FDR"]
    }
       
    ## global difference between stages, resampled reads
    
    pathde=paste(pathDiffExp, sp, "/DifferentialExpression_AllAges_ResampledReads_PCLncRNAs_",tissue,".txt", sep="")

    if(file.exists(pathde)){
      this.de.resampled=read.table(pathde, h=T, stringsAsFactors=F, sep="\t", quote="\"")
      rownames(this.de.resampled)=this.de.resampled$GeneID
      
      de.global.resampled[[paste("FDR", tissue, sep=".")]]=this.de.resampled[allgenes, "FDR"]
    }
    
    ## consecutive stages

    for(i in 1:4){
      stage1=age.order[i]
      stage2=age.order[i+1]
      
      pathde=paste(pathDiffExp, sp, "/DifferentialExpression_ConsecutiveAges_AllReads_PCLncRNAs_", tissue, "_", stage1, "_", stage2, ".txt", sep="")
      
      if(file.exists(pathde)){
        this.de.stages=read.table(pathde, h=T, stringsAsFactors=F, sep="\t", quote="\"")
        rownames(this.de.stages)=this.de.stages$GeneID
      
        de.consecutive.stages[[paste("FDR", paste(tissue, stage1, stage2, sep="_"), sep=".")]]=this.de.stages[allgenes, "FDR"]
      }      
    }
  }

  ## write results for DE global all reads

  if(length(de.global.allreads)>0){
    de.global.allreads=as.data.frame(de.global.allreads, stringsAsFactors=F)
    de.global.allreads$GeneID=allgenes
    de.global.allreads$GeneType=genetype

    de.global.allreads=de.global.allreads[,c("GeneID", "GeneType", setdiff(colnames(de.global.allreads), c("GeneID", "GeneType")))]
    
    write.table(de.global.allreads, file=paste(pathDatasets, "SupplementaryDataset4/DifferentialExpression_GlobalAgeEffect_AllReads_",sp,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  }
  
  ## write results for DE global resampled

  if(length(de.global.resampled)>0){
    de.global.resampled=as.data.frame(de.global.resampled, stringsAsFactors=F)
    de.global.resampled$GeneID=allgenes
    de.global.resampled$GeneType=genetype
    
    de.global.resampled=de.global.resampled[,c("GeneID", "GeneType", setdiff(colnames(de.global.resampled), c("GeneID", "GeneType")))]
    
    write.table(de.global.resampled, file=paste(pathDatasets, "SupplementaryDataset4/DifferentialExpression_GlobalAgeEffect_ResampledReads_",sp,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  }

  ## write results for DE consecutive stages
  
  if(length(de.consecutive.stages)>0){
    de.consecutive.stages=as.data.frame(de.consecutive.stages, stringsAsFactors=F)
    de.consecutive.stages$GeneID=allgenes
    de.consecutive.stages$GeneType=genetype
    
    de.consecutive.stages=de.consecutive.stages[,c("GeneID", "GeneType", setdiff(colnames(de.consecutive.stages), c("GeneID", "GeneType")))]
    
    write.table(de.consecutive.stages, file=paste(pathDatasets, "SupplementaryDataset4/DifferentialExpression_ConsecutiveStages_",sp,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  }
}

########################################################################################


