########################################################################

path="LncEvoDevo/"
pathGeneFamilies=paste(path, "results/ortho_genes/whole_genome_alignments/", sep="")
pathLncRNADataset=paste(path, "results/lncRNA_dataset/", sep="")

release=94
splist=c("Mouse", "Rat", "Chicken")
type=paste("FilteredTranscripts_StringTie_Ensembl",release,sep="")
shorttype="StringTie"

########################################################################

genefamilies=read.table(paste(pathGeneFamilies, "GeneFamilies_1to1_",shorttype,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

########################################################################

for(sp in splist){
  allinfo=read.table(paste(pathLncRNADataset, sp, "/AllInfo_WithRejectionReason_",type,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(allinfo)=allinfo$GeneID

  gt=rep("other", dim(allinfo)[1])
  names(gt)=rownames(allinfo)

  gt[which(allinfo$EnsemblBiotype=="protein_coding")]="protein_coding"
  gt[which(allinfo$RejectionReason=="")]="lncRNA"

  rr=allinfo$RejectionReason
  names(rr)=rownames(allinfo)

  rr[which(rr=="")]="none"

  genefamilies[,paste("GeneType.",sp,sep="")]=rep(NA, dim(genefamilies)[1])
  genefamilies[which(!is.na(genefamilies[,paste("ID.",sp,sep="")])),paste("GeneType.",sp,sep="")]=gt[genefamilies[which(!is.na(genefamilies[,paste("ID.",sp,sep="")])),paste("ID.",sp,sep="")]]
  
  genefamilies[,paste("RejectionReason.",sp,sep="")]=rep(NA, dim(genefamilies)[1])
  genefamilies[which(!is.na(genefamilies[,paste("ID.",sp,sep="")])),paste("RejectionReason.",sp,sep="")]=rr[genefamilies[which(!is.na(genefamilies[,paste("ID.",sp,sep="")])),paste("ID.",sp,sep="")]]
  
}

write.table(genefamilies, file=paste(pathLncRNADataset,"/GeneFamilies_",type,".txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

########################################################################


