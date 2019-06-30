##########################################################################

path="LncEvoDevo/"
pathResults=paste(path, "results/intergenic_regions/", sep="")

release=94

options(scipen=999)

##########################################################################

for(sp in c("Mouse", "Rat")){
  ig=read.table(paste(pathResults, sp, "/ResampledIntergenicRegions_Size5kb_MinDistance5kb_FilteredTranscripts_StringTie_Ensembl94.txt", sep=""), h=F, stringsAsFactors=F)

  ig[,2]=format(ig[,2], trim=TRUE, scientific=FALSE)
  ig[,3]=format(ig[,3], trim=TRUE, scientific=FALSE)

  ig=ig[which(ig[,1]%in%c(as.character(1:19), "X", "Y")),]

  ig$ID=apply(ig, 1, function(x) paste(x, collapse=","))
  ig[,1]=paste("chr", ig[,1], sep="")

  write.table(ig, file=paste(pathResults, sp, "/ResampledIntergenicRegions_Size5kb_MinDistance5kb_FilteredTranscripts_StringTie_Ensembl94.bed", sep=""), row.names=F, col.names=F, sep="\t", quote=F)
}

##########################################################################

