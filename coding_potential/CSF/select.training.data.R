##########################################################################

path="LncEvoDevo/"
pathOrtho=paste(path,"data/ensembl_ortho/",sep="")
pathCSF=paste(path,"results/CSF/",sep="")
pathAnnot=paste(path,"data/ensembl_annotations/",sep="")

##########################################################################
#c("Mouse", "Rat", "Chicken")
for(sp in c("Human")){
 
  print(sp)
  
  ##########################################################################
  
  gene.info=read.table(paste(pathAnnot, sp,"/GeneInfo_Ensembl84.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  tx.info=read.table(paste(pathAnnot, sp,"/TranscriptInfo_Ensembl84.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  
  gene.info=gene.info[which(gene.info$biotype=="protein_coding" & gene.info$status=="KNOWN" & gene.info$name%in%c(as.character(1:30), "X", "Y", "Z", "W")),]
  tx.info=tx.info[which(tx.info$biotype=="protein_coding" & tx.info[,1]%in%gene.info[,1]),]
  
  ## select longest transcript (in terms of genomic coordinates)
  
  tx.info$Length=tx.info$seq_region_end-tx.info$seq_region_start+1
  
  max.tx.length=tapply(tx.info$Length, as.factor(tx.info[,1]), max)
  names(max.tx.length)=levels(as.factor(tx.info[,1]))
  tx.info$MaxLength=max.tx.length[tx.info[,1]]
  
  tx.info=tx.info[which(tx.info$Length==tx.info$MaxLength),]
  dupli=which(duplicated(tx.info[,1]))
  
  if(length(dupli)>0){
    tx.info=tx.info[-dupli,]
  }
  
  ##########################################################################
  
  ## select 3,000 genes with 1-to-1 ortho in mouse, rat, human, opossum, chicken
  
  ortho=read.table(paste(pathOrtho, "GeneFamilies_1to1_AllSpecies_ProteinCoding_Ensembl84.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

  ortho=ortho[which((!is.na(ortho$Mouse)) & (!is.na(ortho$Rat)) & (!is.na(ortho$Chicken)) & (!is.na(ortho$Human)) & (!is.na(ortho$Opossum)) & ortho[,sp]%in%tx.info[,1]),]
  
  
  selected.genes=ortho[,sp][1:3000]
  
  tx.info=tx.info[which(tx.info[,1]%in%selected.genes),]
  selected.transcripts=tx.info[,2]
  
  writeLines(selected.genes, con=paste(pathCSF, "training_dataset/",sp,"/SelectedGenes.txt",sep=""))
  writeLines(selected.transcripts, con=paste(pathCSF, "training_dataset/",sp,"/SelectedTranscripts.txt",sep=""))
  
 ##########################################################################

  introns=read.table(paste(pathCSF,"intron_regions/",sp,"/IntronRegions_MinSize500_MaxSize10000_ExcludeLength200_Ensembl84.txt",sep=""), h=F, stringsAsFactors=F, sep="\t")
  
  introns=introns[which(introns$V1%in%selected.genes),]
  
  introns$V3=paste("chr", introns$V3, sep="")
  
  results.int=introns[,c(1,3,4,5,6)]
  results.int[,6]=rep(0,dim(results.int)[1])
  
  colnames(results.int)=c("Gene", "Chr", "Start", "End", "Strand", "Frame")
  
  write.table(results.int, file=paste(pathCSF, "training_dataset/",sp,"/SelectedIntrons.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  
  ##########################################################################
  
  cds.coords=read.table(paste(pathAnnot, sp,"/CDSCoords_Ensembl84.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")

  cds.coords=cds.coords[which(cds.coords$Ensembl.Gene.ID%in%selected.genes & cds.coords$Ensembl.Transcript.ID%in%selected.transcripts & cds.coords$Exon.Chr.Start..bp.==cds.coords$Genomic.coding.start & cds.coords$Exon.Chr.End..bp.==cds.coords$Genomic.coding.end),]
  
  results.cds=cds.coords[,c("Ensembl.Gene.ID", "Chromosome.Name", "Genomic.coding.start", "Genomic.coding.end", "Strand", "start.phase")]
  
  colnames(results.cds)=c("Gene", "Chr", "Start", "End", "Strand", "Frame")
  
  results.cds[which(results.cds$Strand==-1),"Frame"]=cds.coords[which(results.cds$Strand==-1),"end.phase"]

  nbrev=length(which(results.cds$Strand==-1))

  print(paste("changed strand for ",nbrev, "exons"))
  
  results.cds$Chr=paste("chr",results.cds$Chr,sep="")
  
  write.table(results.cds, file=paste(pathCSF, "training_dataset/",sp,"/SelectedCDS.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  
##########################################################################
  
}
