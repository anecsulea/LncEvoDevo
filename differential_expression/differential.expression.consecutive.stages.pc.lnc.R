########################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathResults=paste(path, "results/differential_expression/", sep="")

set.seed(19)

ages.order=c("EarlyEmbryo", "LateEmbryo", "Newborn", "Adult", "Aged")

library(DESeq2)

## R3.5.0
## DESeq2_1.22.1   

splist=c("Mouse", "Rat", "Chicken")

########################################################################

for(sp in splist){
  print(sp)

  load(paste("RData/reads.resampling.pc.lncRNA.",sp,".RData",sep=""))

  reads=rbind(reads.pc, reads.lnc)
     
  for(tiss in c("Brain", "Kidney", "Liver", "Testis")){
    
    if(tiss=="Testis"){
      firststage=2
    } else{
      firststage=1
    }
    
    for(i in firststage:(length(ages.order)-1)){
      age1=ages.order[i]
      age2=ages.order[i+1]

      if(length(which(tissues==tiss & ages==age1))>0 & length(which(tissues==tiss & ages==age2))>0){
        print(paste(tiss, age1, age2))
        
        this.reads=reads[,which(tissues==tiss & ages%in%c(age1, age2))]
        this.ages=as.factor(ages[which(tissues==tiss & ages%in%c(age1, age2))])
        
        colnames(this.reads)=paste("NbReads", colnames(this.reads), sep=".")
        
        colData=data.frame("condition"=this.ages)
        dds=DESeqDataSetFromMatrix(countData = this.reads, colData=colData, design = ~ condition)
        dds=DESeq(dds,test="Wald")
        
        coef1=paste("condition_",age2,"_vs_",age1,sep="")
        coef2=paste("condition_",age1,"_vs_",age2,sep="")
        sign=1
        
        if(coef1%in%resultsNames(dds)){
          coef=coef1
        } else{
          if(coef2%in%resultsNames(dds)){
            coef=coef2
          sign=-1
          } else{
            stop("cannot find coefficient!")
          }
        }
        
        res=lfcShrink(dds, coef=coef, type="apeglm") ## shrunk log2 fold change
        
        res=res[order(res$pvalue),]
        res=as.data.frame(res)
        
        res=res[,c("pvalue", "padj", "log2FoldChange")]
        colnames(res)=c("PValue", "FDR", "log2FoldChange")
        
        res$log2FoldChange=sign * res$log2FoldChange ## always age2 vs age1
        
        res$GeneID=rownames(res)
        res=cbind(res, this.reads[rownames(res),])
        
        res=res[,c("GeneID", setdiff(colnames(res), "GeneID"))]
        
        nbsignif1=length(which(res$FDR<0.1))
        nbsignif2=length(which(res$FDR<0.01))
        
        print(paste(nbsignif1, "significant FDR < 0.1"))
        print(paste(nbsignif2, "significant FDR < 0.01"))
        
        write.table(res, file=paste(pathResults, sp, "/DifferentialExpression_ConsecutiveAges_AllReads_PCLncRNAs_",tiss, "_",age1,"_",age2,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
      }
    }
  }
}

########################################################################
