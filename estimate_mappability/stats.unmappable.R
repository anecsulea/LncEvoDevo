###############################################################

path="/pandata/necsulea/LncEvoDevo/"
pathMap=paste(path, "results/mappability/", sep="")

###############################################################

for(sp in c("Mouse", "Rat", "Chicken")){
  print(sp)
  
  unmap=read.table(paste(pathMap, sp, "/unmappable_regions.txt", sep=""), h=F, stringsAsFactors=F)
  chrsizes=read.table(paste(pathMap, sp, "/chromosome_sizes_samtools.txt", sep=""), h=F, stringsAsFactors=F)

  chrsizes=chrsizes[which(chrsizes$V1%in%unique(unmap$V1)),]
  totsize=sum(as.numeric(chrsizes$V2))

  totunmap=sum(as.numeric(unmap$V3)-as.numeric(unmap$V2)+1)

  pcunmap=round(100*totunmap/totsize, digits=2)
  totmap=totsize-totunmap
  
  print(paste(totmap,"mappable bases"))
  print(paste(pcunmap,"% unmappable bases"))
  
}

###############################################################


