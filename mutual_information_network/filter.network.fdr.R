############################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathResults=paste(path, "results/mutual_information_network/", sep="")

options(scipen=999)

maxFDR=0.001

############################################################################

for(type in c("lncRNAs_only")){ ##, "all_genes"
  for(sp in c("Mouse", "Rat")){
    
    network=read.table(paste(pathResults, type, "/", sp,"/network.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    network$FDR=p.adjust(network$pvalue, method="BH")
    
    filtered=network[which(network$FDR<maxFDR),]

    write.table(filtered, file=paste(pathResults, type, "/", sp,"/filtered_network_FDR",maxFDR,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
 
  }
}

############################################################################
