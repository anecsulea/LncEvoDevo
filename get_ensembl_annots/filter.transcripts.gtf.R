##############################################################################

path="LncEvoDevo/"
pathAnnot=paste(path, "data/ensembl_annotations/",sep="")
pathIndexes=paste(path, "data/genome_indexes/",sep="")

release=94

options(scipen=999) ## remove scientific notation ## options(scipen=0) to get it back

maxlen=2500000

##############################################################################

for(sp in c("Human", "Mouse", "Rat", "Chicken")){ 

  print(sp)

  #################################

  ## readthrough transcripts

  rt=read.table(paste(pathAnnot, sp, "/ReadthroughTranscripts_Ensembl",release, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

  #################################

  txinfo=read.table(paste(pathAnnot, sp, "/TranscriptInfo_Ensembl",release,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  colnames(txinfo)[1]="gene"
  colnames(txinfo)[2]="tx"

  txinfo$Length=txinfo[,"seq_region_end"]-txinfo[,"seq_region_start"]+1
  oklen=txinfo$tx[which(txinfo$Length<=maxlen)]
  notok=setdiff(txinfo$tx, oklen)

  if(length(notok)>0){
    print(paste("Discarded long transcripts:",paste(notok, collapse=", "))) 
  }

  geneinfo=read.table(paste(pathAnnot, sp, "/GeneInfo_Ensembl",release,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

  pcgenes=geneinfo$stable_id[which(geneinfo$biotype=="protein_coding")]
  pctx=txinfo[which(txinfo$biotype=="protein_coding"),2]

  #################################

  ## for human, some genes have transcripts on multiple strands, we remove them

  nbstrands=tapply(txinfo$seq_region_strand, as.factor(txinfo$gene), function(x) length(unique(x)))
  weirdgenes=names(nbstrands)[which(nbstrands!=1)]

  if(length(weirdgenes)>0){
    print(paste(length(weirdgenes), "genes with more than one strand"))
    print(weirdgenes)
  }

  #################################

  print(paste(length(which(txinfo$tx%in%rt$TranscriptID)), " read-through transcripts"))
  
  ok=which(((txinfo$gene%in%pcgenes & txinfo$tx%in%pctx) | (!(txinfo$gene%in%pcgenes))) & txinfo$tx%in%oklen & (!(txinfo$tx%in%rt$TranscriptID)) & (!txinfo$gene%in%weirdgenes))
 
  selected=txinfo[ok,2]

  print(paste(length(selected), "selected transcripts"))

  #################################

  print("reading gtf")
  
  gtf=read.table(paste(pathAnnot, sp, "/AllTranscripts_Ensembl",release,".gtf",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="")
  
  colnames(gtf)=c("chr", "source", "type", "start", "end", "score1", "strand", "score2", "info")
  
  gtf=gtf[which(gtf$type=="exon"),]

  print("done")
  
  ##################################
  
  gtf.info=lapply(gtf$info, function(x) unlist(strsplit(x, split=";")))
  txid=unlist(lapply(gtf.info, function(x) grep("transcript_id", x, value=T)))
  txid=unlist(lapply(txid, function(x) unlist(strsplit(x,split="\""))[2]))
  
  ###################################

  gtf.selected=gtf[which(txid%in%selected),]
  
  write.table(gtf.selected, file=paste(pathAnnot, sp,"/FilteredTranscripts_Ensembl",release,".gtf",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

  ###################################
}


##############################################################################
