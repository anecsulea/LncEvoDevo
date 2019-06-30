############################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathResults=paste(path, "results/mutual_information_network/", sep="")

options(scipen=999)

maxFDR=0.001

############################################################################

type="lncRNAs_only"

network=read.table(paste(pathResults, type, "/common_interactions_FDR",maxFDR,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")


############################################################################

for(sp in c("Mouse", "Rat")){

  print(sp)
 
  tpm=read.table(paste(pathResults, type, "/",sp,"/TPM.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

  rownames(tpm)=tpm$Gene
  tpm=tpm[,which(colnames(tpm)!="Gene")]

  tpm=log2(tpm+1)
  tpm=as.matrix(tpm)

  cormat.pearson=cor(t(tpm), method="pearson")
  cormat.spearman=cor(t(tpm), method="spearman")

  print("pearson")
  corPearson=unlist(lapply(1:dim(network)[1], function(x) cormat.pearson[network$ID1[x], network$ID2[x]]))

  print("spearman")
  corSpearman=unlist(lapply(1:dim(network)[1], function(x) cormat.spearman[network$ID1[x], network$ID2[x]]))
  
  network[, paste("PearsonCorr.",sp,sep="")]=corPearson
  network[, paste("SpearmanCorr.",sp,sep="")]=corSpearman
}

############################################################################

write.table(network, file=paste(pathResults, type, "/common_interactions_with_correlation_values_FDR",maxFDR,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

############################################################################

for(method in c("Pearson", "Spearman")){
  pdf(file=paste("Correlation_vs_MutualInformation_",method,".pdf", sep=""), width=12, height=6)

  par(mfrow=c(1,2))

  smoothScatter(network[,paste(method, "Corr.Mouse",sep="")], network$MI.Mouse, pch=20, xlab=paste(method, "'s correlation",sep=""), ylab="Mutual information",  main="mouse")

  smoothScatter(network[,paste(method, "Corr.Rat",sep="")], network$MI.Rat, pch=20, xlab=paste(method, "'s correlation",sep=""), ylab="Mutual information",  main="rat")
  
  dev.off()
}

############################################################################

for(method in c("Pearson", "Spearman")){
  pdf(file=paste("Correlation_",method,"_SignificantInteractions.pdf", sep=""), width=12, height=6)

  par(mfrow=c(1,2))

  plot(density(network[,paste(method, "Corr.Mouse",sep="")]),  main="mouse", xlab=paste(method, "'s correlation",sep=""), ylab="Density")
  
  
  plot(density(network[,paste(method, "Corr.Rat",sep="")]),  main="rat", xlab=paste(method, "'s correlation",sep=""), ylab="Density")

  dev.off()
}

############################################################################
