##############################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/",sep="")
pathTables=paste(path, "supplementary_tables/",sep="")
pathMainFigures=paste(path, "scripts/main_figures/",sep="")

##############################################################################

source(paste(pathMainFigures, "compute.go.enrichment.R",sep=""))

##############################################################################

tissue.order=c("Brain","Kidney","Liver","Testis")
age.order=c("EarlyEmbryo", "LateEmbryo", "Newborn", "Adult", "Aged")
newage=c("MidStageEmbryo","LateEmbryo", "Newborn", "Adult", "Aged")
names(newage)=age.order

maxFDR=0.01

##############################################################################

for(sp in c("Mouse", "Rat")){
  load(paste(pathMainFigures,"RData/data.gene.ontology.",sp,".RData", sep=""))
  
  diffexp=read.table(paste(pathDatasets, "SupplementaryDataset4/DifferentialExpression_ConsecutiveStages_",sp,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  rownames(diffexp)=diffexp$GeneID
  diffexp=diffexp[which(diffexp$GeneType=="ProteinCoding"),]
  
  tpm=read.table(paste(pathDatasets, "SupplementaryDataset3/Statistics_AverageTPM_",sp,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
  rownames(tpm)=tpm$GeneID
  tpm=tpm[rownames(diffexp),]

  for(tiss in tissue.order){
    this.age.order=age.order
    if(tiss=="Testis"){
      this.age.order=setdiff(this.age.order, "EarlyEmbryo")
    }
    
    for(i in 1:(length(this.age.order)-1)){
      stage1=this.age.order[i]
      stage2=this.age.order[i+1]

      print(paste(sp, tiss, stage1, stage2))
      
      up=diffexp$GeneID[which(diffexp[,paste("FDR.",tiss,"_",stage1,"_",stage2,sep="")]<maxFDR & tpm[,paste("MeanTPM.",tiss, "_",stage2,sep="")]>tpm[,paste("MeanTPM.",tiss, "_",stage1,sep="")])]
      
      down=diffexp$GeneID[which(diffexp[,paste("FDR.",tiss,"_",stage1,"_",stage2,sep="")]<maxFDR & tpm[,paste("MeanTPM.",tiss, "_",stage2,sep="")]<tpm[,paste("MeanTPM.",tiss, "_",stage1,sep="")])]
      
      all=diffexp$GeneID[which(!is.na(diffexp[,paste("FDR.",tiss,"_",stage1,"_",stage2,sep="")]))]
      
      go.up=compute.go.enrichment(up, all, GOdata[["biological_process"]][["categories"]], GOdata[["biological_process"]][["genelist"]])
      go.up=go.up[order(go.up$FDR),]
      go.up=go.up[which(go.up$FDR<0.1),]
      
      if(dim(go.up)[1]>0){
        write.table(go.up, file=paste(pathDatasets, "SupplementaryDataset4/GOEnrichment_DifferentialExpression_ConsecutiveStages/BiologicalProcess_UpregulatedGenes_",sp,"_",tiss, "_",newage[stage1],"_",newage[stage2],".txt",sep=""),row.names=F, col.names=T, sep="\t", quote=F)
      }
      
      go.down=compute.go.enrichment(down, all, GOdata[["biological_process"]][["categories"]], GOdata[["biological_process"]][["genelist"]])
      go.down=go.down[order(go.down$FDR),]
      go.down=go.down[which(go.down$FDR<0.1),]
      
      if(dim(go.down)[1]>0){
        write.table(go.down, file=paste(pathDatasets, "SupplementaryDataset4/GOEnrichment_DifferentialExpression_ConsecutiveStages/BiologicalProcess_DownregulatedGenes_",sp,"_",tiss, "_",newage[stage1],"_",newage[stage2],".txt",sep=""),row.names=F, col.names=T, sep="\t", quote=F)
      }        
      
    }
  }

}

##############################################################################
