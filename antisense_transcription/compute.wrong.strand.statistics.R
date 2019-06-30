#########################################################################

path="LncEvoDevo/"
pathHisat=paste(path, "results/hisat/", sep="")
pathEnsembl=paste(path, "data/ensembl_annotations/", sep="")

splist=c("Mouse", "Rat", "Chicken")

options(scipen=999) ## remove scientific notation ## options(scipen=0) to get it back

#########################################################################

for(sp in splist){

  print(sp)

  annotated=read.table(paste(pathEnsembl, sp, "/SpliceSites_Ensembl89.txt", sep=""), h=F, stringsAsFactors=F)
  annotated$strand=rep(NA, dim(annotated)[1])
  annotated$strand[which(annotated$V4=="+")]="1"
  annotated$strand[which(annotated$V4=="-")]="-1"

  colnames(annotated)=c("chr", "start", "end", "oldstrand", "strand")
  
  id.annotated=paste(annotated$chr, annotated$start+2, annotated$end, annotated$strand, sep=",")


  samples=system(paste("ls ", pathHisat, sp, "/ | grep -v txt", sep=""), intern=T)

  stats=list()

  for(sample in samples){
    print(sample)
    
    correct=read.table(paste(pathHisat, sp, "/", sample, "/junctions.txt", sep=""), h=T, stringsAsFactors=F)
    wrong=read.table(paste(pathHisat, sp, "/", sample, "/junctions_wrongstrand_counts.txt", sep=""), h=T, stringsAsFactors=F)

    correct$id=paste(correct$Chromosome, correct$Start, correct$End, correct$ProbableStrand, sep=",")
    
    wrong$Start=unlist(lapply(wrong$Coords, function(x) unlist(strsplit(x, split="_"))[1]))
    wrong$End=unlist(lapply(wrong$Coords, function(x) unlist(strsplit(x, split="_"))[2]))
    
    wrong$id=paste(wrong$Chr, wrong$Start, wrong$End, wrong$SpliceStrand, sep=",")

    correct=correct[which(correct$id%in%id.annotated),]
    wrong=wrong[which(wrong$id%in%id.annotated),]
    
    nbreads.correct=sum(correct$NbReads)
    nbreads.wrong=sum(wrong$NbReads)

    stats[[sample]]=c(nbreads.correct, nbreads.wrong)
  }

  stats=t(as.matrix(as.data.frame(stats)))
  stats=as.data.frame(stats)
  colnames(stats)=c("NbCorrect", "NbWrong")
  rownames(stats)=samples

  stats$PropWrong=stats$NbWrong/(stats$NbCorrect+stats$NbWrong)
  stats$SampleID=rownames(stats)
  stats=stats[,c("SampleID", "NbCorrect", "NbWrong", "PropWrong")]
  write.table(stats, file=paste(pathHisat, sp, "/stats_antisense_junctions_annotated.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

  
}

#########################################################################


