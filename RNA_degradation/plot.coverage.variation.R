################################################################################

path="LncEvoDevo/"
pathCoverage=paste(path,"results/RNA_degradation/",sep="")
pathEnsembl=paste(path,"data/ensembl_annotations/",sep="")

################################################################################

for(sp in c("Mouse", "Rat", "Chicken")){
  samples=system(paste("ls ",pathCoverage, sp, "/ | grep -v txt", sep=""),intern=T)
  print(samples)
  
  ################################################################################
  
  geneinfo=read.table(paste(pathEnsembl,sp, "/GeneInfo_Ensembl94.txt",sep=""),h=F,skip=1,sep="\t",stringsAsFactors=F,quote="\"")
  colnames(geneinfo)=c("GeneID","Biotype","Description","Chr","Start","End","Strand")
  rownames(geneinfo)=geneinfo$GeneID
  
  pc=geneinfo$GeneID[which(geneinfo$Biotype=="protein_coding")]
  
################################################################################
  
  output=file(paste(pathCoverage,sp, "/degradation_stats.txt",sep=""),open="w")

  writeLines(paste(c("SampleID","Slope", "PValue","CoverageBias"), collapse="\t"), output)
  
  for(sample in samples){
    print(paste(sp, sample))
    
    if(file.exists(paste(pathCoverage,sp, "/", sample,"/CoverageVariation_Ensembl94_20windows_minlength400.txt",sep=""))){
      cov=read.table(paste(pathCoverage,sp, "/", sample,"/CoverageVariation_Ensembl94_20windows_minlength400.txt",sep=""),h=T,stringsAsFactors=F)
      covpc=cov[which(cov$GeneID%in%pc),6:25]
      covpc=as.matrix(covpc)
      
      meancov=apply(covpc,2,mean,na.rm=T)
      mediancov=apply(covpc,2,median,na.rm=T)
      
      pos=5:17
      lm1=lm(meancov[5:17]~pos)
      
      slope=summary(lm1)$coefficients[2,1]
      pval=summary(lm1)$coefficients[2,4]

      ratio=(meancov[17]-meancov[5])/meancov[17]
      
      writeLines(paste(sample, slope, pval, ratio, sep="\t"),output)
      
      pdf(file=paste("figures/",sp, "_", sample,"_MedianCoverage.pdf",sep=""))
      par(mar=c(5.1,5.1,4.1,2.1))
      plot(1:20,mediancov,pch=20,type="b",axes=F,xlab="Position along gene length",ylab="Median read coverage",cex.lab=1.5,main=sample,cex.main=1.5)
      axis(side=1,at=seq(from=2,to=20,by=2),labels=seq(from=10,to=100,by=10),cex.axis=1.25)
      axis(side=2,cex.axis=1.25)
      box()
      dev.off()
      
      pdf(file=paste("figures/",sp, "_", sample,"_MeanCoverage.pdf",sep=""))
      par(mar=c(5.1,5.1,4.1,2.1))
      plot(1:20,meancov,pch=20,type="b",axes=F,xlab="Position along gene length",ylab="Mean read coverage",cex.lab=1.5,main=sample,cex.main=1.5)
      axis(side=1,at=seq(from=2,to=20,by=2),labels=seq(from=10,to=100,by=10),cex.axis=1.25)
      axis(side=2,cex.axis=1.25)
      box()
      dev.off()
      
    }
  }
  
  close(output)
}
################################################################################
