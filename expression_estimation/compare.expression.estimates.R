#################################################################################

path="LncEvoDevo/"
pathExpression=paste(path,"results/expression_estimation/", sep="")
pathAllInfo=paste(path,"results/lncRNA_dataset/", sep="")

types=c("StringTie")
splist=c("Mouse", "Rat", "Chicken")

#################################################################################

for(sp in splist){
  allinfo=read.table(paste(pathAllInfo, sp, "/AllInfo_WithRejectionReason_FilteredTranscripts_StringTie_Ensembl94.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  pc=allinfo$GeneID[which(allinfo[,"EnsemblBiotype"]=="protein_coding")]
  lnc=allinfo$GeneID[which(allinfo[,"RejectionReason"]=="")]
  
  for(type in types){
    print(paste(sp, type))
    
    subread=read.table(paste(pathExpression, sp, "/AllSamples_UniqueReadCounts_",type,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    rownames(subread)=subread$GeneID
    nbgenes.subread=dim(subread)[1]

    print(paste(nbgenes.subread, "genes with subread"))
    
    kallisto=read.table(paste(pathExpression, sp, "/AllSamples_KallistoEstimatedCounts_",type,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    rownames(kallisto)=kallisto$GeneID
    nbgenes.kallisto=dim(kallisto)[1]

    print(paste(nbgenes.kallisto, "genes with Kallisto"))

    genes=intersect(rownames(kallisto), rownames(subread))
    subread=subread[genes,]
    kallisto=kallisto[genes,]

    pc=intersect(pc, genes)
    lnc=intersect(lnc, genes)
    
    samples=setdiff(colnames(kallisto), "GeneID")

    for(sample in samples){

      ## all genes
      
      R=cor(subread[,sample], kallisto[,sample], method="spearman")

      png(file=paste("figures/ExpressionCorrelations_",sp,"_",type,"_",sample,"_AllGenes.png",sep=""))
      
      plot(log2(subread[,sample]+1), log2(kallisto[,sample]+1), xlab="subread", ylab="kallisto", main=paste(sp, sample, "R =", round(R, digits=2)), pch=20)

      dev.off()

      ## pc genes only

      R=cor(subread[pc,sample], kallisto[pc,sample], method="spearman")

      png(file=paste("figures/ExpressionCorrelations_",sp,"_",type,"_",sample,"_ProteinCoding.png",sep=""))
      
      plot(log2(subread[pc,sample]+1), log2(kallisto[pc,sample]+1), xlab="subread", ylab="kallisto", main=paste(sp, sample, "R =", round(R, digits=2)), pch=20)

      dev.off()

      ## lnc genes only

      R=cor(subread[lnc,sample], kallisto[lnc,sample], method="spearman")

      png(file=paste("figures/ExpressionCorrelations_",sp,"_",type,"_",sample,"_LncRNAs.png",sep=""))
      
      plot(log2(subread[lnc,sample]+1), log2(kallisto[lnc,sample]+1), xlab="subread", ylab="kallisto", main=paste(sp, sample, "R =", round(R, digits=2)), pch=20)

      dev.off()
    }
  }
}

#################################################################################

