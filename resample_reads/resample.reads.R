############################################################################

path="LncEvoDevo/"
pathHisat=paste(path, "results/hisat/", sep="")

set.seed(19)
options(scipen=999)

##########################################################################

clean=TRUE

args <- commandArgs(trailingOnly = TRUE)

if(length(args)>0){
  for(k in 1:length(args)){
    eval(parse(text=args[[k]]))
  }
}

print(paste("sp=",sp))
print(paste("sample=",sample))

############################################################################

nb.resampled=read.table(paste(pathHisat, sp, "/nb_resampled_unique_reads_noMT.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

if(sample%in%nb.resampled$Sample){
  nbreads=nb.resampled$NbReads[which(nb.resampled$Sample==sample)]
  resampled=nb.resampled$NbResampled[which(nb.resampled$Sample==sample)]
  
  print(paste(nbreads, resampled))
  
  
  if(file.exists(paste(pathHisat, sp, "/", sample, "/resampled_read_ids_noMT/read_ids_", resampled,".txt",sep=""))){
    print(paste("already done",resampled))
  } else{
    this.sample=sample(0:(nbreads-1), size=resampled, replace=F)
    writeLines(as.character(this.sample), con=paste(pathHisat, sp, "/", sample, "/resampled_read_ids_noMT/read_ids_", resampled,".txt",sep=""))
  }
} else{
  print(paste("not doing resampling for ",sp, sample))
}
  
############################################################################
