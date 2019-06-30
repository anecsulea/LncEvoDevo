#########################################################################

path="LncEvoDevo/"
pathGeneOverlaps=paste(path, "results/gene_overlaps/", sep="")

release=94

#########################################################################

for(sp in c("Mouse", "Rat", "Chicken")){
  exons=read.table(paste(pathGeneOverlaps, sp, "/ExonBlocks_ExcludingOverlapOtherGenes_FilteredTranscripts_StringTie_Ensembl",release,".txt",sep=""), h=F, stringsAsFactors=F)

  colnames(exons)=c("GeneID", "ExonID", "Chr", "Start", "End", "Strand")

  exons$Strand=as.character(exons$Strand)
  exons$Strand[which(exons$Strand=="1")]="+"
  exons$Strand[which(exons$Strand=="-1")]="-"
                    

  info=paste("gene_id \"",exons$GeneID, "\";", sep="")
  
  writeLines(paste(exons$Chr, rep("ExonsNoOverlap", dim(exons)[1]), rep("exon", dim(exons)[1]), exons$Start, exons$End, rep(".", dim(exons)[1]), exons$Strand, rep(".", dim(exons)[1]), info, sep="\t"), con=paste(pathGeneOverlaps, sp, "/ExonBlocks_ExcludingOverlapOtherGenes_FilteredTranscripts_StringTie_Ensembl",release,".gtf",sep=""))
}

#########################################################################


