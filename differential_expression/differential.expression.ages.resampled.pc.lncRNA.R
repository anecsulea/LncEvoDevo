########################################################################

path="LncEvoDevo/"
pathResults=paste(path, "results/differential_expression/", sep="")

########################################################################

set.seed(19)

options(stringsAsFactors=F)

library(DESeq2)

########################################################################

for(sp in c("Mouse", "Rat")){
  print(sp)

  load(paste("RData/reads.resampling.pc.lncRNA.",sp,".RData",sep=""))
  
  print("resampled")

  reads=rbind(reads.pc.resampled, reads.lnc)

  for(tiss in c("Brain", "Kidney", "Liver", "Testis")){
    print(tiss)
    
    this.reads=reads[,which(tissues==tiss)]
    this.ages=as.factor(ages[which(tissues==tiss)])

    print(this.ages)
    
    colnames(this.reads)=paste("NbReads", colnames(this.reads), sep=".")
      
    colData=data.frame("condition"=this.ages)
    dds=DESeqDataSetFromMatrix(countData = this.reads, colData=colData, design = ~ condition)
    dds=DESeq(dds,test="LRT", reduced = ~ 1, minReplicatesForReplace=15 )

    res=results(dds)
    res=res[order(res$pvalue),]
    res=as.data.frame(res)

    res=res[,c("pvalue", "padj")]
    colnames(res)=c("PValue", "FDR")
    
    res$GeneID=rownames(res)
    res=cbind(res, this.reads[rownames(res),])
    
    res=res[,c("GeneID", setdiff(colnames(res), "GeneID"))]

    write.table(res, file=paste(pathResults, sp, "/DifferentialExpression_AllAges_ResampledReads_PCLncRNAs_",tiss, ".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  }

  print("not resampled")
  
  reads=rbind(reads.pc, reads.lnc)

  for(tiss in c("Brain", "Kidney", "Liver", "Testis")){
    print(tiss)
    
    this.reads=reads[,which(tissues==tiss)]
    this.ages=as.factor(ages[which(tissues==tiss)])
    
    colnames(this.reads)=paste("NbReads", colnames(this.reads), sep=".")
      
    colData=data.frame("condition"=this.ages)
    dds=DESeqDataSetFromMatrix(countData = this.reads, colData=colData, design = ~ condition)
    dds=DESeq(dds,test="LRT", reduced = ~ 1, minReplicatesForReplace=15 )

    res=results(dds)
    res=res[order(res$pvalue),]
    res=as.data.frame(res)

    res=res[,c("pvalue", "padj")]
    colnames(res)=c("PValue", "FDR")
    
    res$GeneID=rownames(res)
    res=cbind(res, this.reads[rownames(res),])
   
    res=res[,c("GeneID", setdiff(colnames(res), "GeneID"))]
        
    write.table(res, file=paste(pathResults, sp, "/DifferentialExpression_AllAges_AllReads_PCLncRNAs_", tiss, ".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  }
}

########################################################################
