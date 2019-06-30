######################################################################

path="LncEvoDevo/"
pathStringTie=paste(path,"results/stringtie_assembly/",sep="")
pathEnsembl=paste(path, "data/ensembl_annotations/",sep="")
pathUCSC=paste(path, "data/UCSC_sequences/",sep="")

release=94

options(scipen=999) ## remove scientific notation ## options(scipen=0) to get it back

######################################################################

for(sp in c("Mouse", "Rat", "Chicken")){

  corresp=read.table(paste(pathUCSC,sp,"/chromosomes_Ensembl_UCSC.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
  rownames(corresp)=corresp[,1]

  ###### only StringTie assembly
  
  print(paste(sp, "StringTie"))
    
  st=read.table(paste(pathStringTie,sp,"/combined/ExonBlocks_FilteredTranscripts_StringTie_Ensembl", release,".txt", sep=""), h=F, stringsAsFactors=F, sep="\t", quote="")

  st=st[,c(3,4,5,6)]
  colnames(st)=c("chr", "start", "end", "strand")

  st$id=paste(st$chr, st$start, st$end, st$strand, sep=",")
 
  print(dim(st)[1])
 
  print(all(st$chr%in%rownames(corresp)))
  st$chr=corresp[st$chr,2]
  print(dim(st)[1])

  st$score=rep("1000", dim(st)[1])
  st$start=st$start-1 ## bed format

  dupli=which(duplicated(st$id))

  if(length(dupli)>0){
    st=st[-dupli,]
    print(paste("removed", length(dupli), "duplicated lines"))
  }

  st$outstrand=rep(NA, dim(st)[1])
  st$outstrand[which(st$strand=="1")]="+"
  st$outstrand[which(st$strand=="-1")]="-"

  print(length(which(is.na(st$outstrand))))
  
  write.table(st[,c("chr", "start", "end", "id", "score", "outstrand")], file=paste(pathStringTie,sp,"/combined/ExonBlocks_FilteredTranscripts_StringTie_Ensembl",release,".bed", sep=""), row.names=F, col.names=F, sep="\t", quote=F)
}

######################################################################

