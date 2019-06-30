########################################################################

path="LncEvoDevo/"
pathExpression=paste(path,"results/read_counts_resampled_noMT/", sep="")
pathAnnot=paste(path,"data/ensembl_annotations/", sep="")
pathHisat=paste(path,"results/hisat/", sep="")
pathStringTie=paste(path,"results/stringtie_assembly/",sep="")

set.seed(19)

options(scipen=999)

########################################################################

types=c("StringTie", "Ensembl")

########################################################################

for(sp in c("Mouse", "Rat", "Chicken")){ 
  nbreads=read.table(paste(pathHisat, sp, "/nb_resampled_unique_reads_noMT.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  samples=nbreads$Sample
  
  for(type in types){

    unique.reads=list()
    
    for(sample in samples){
      nb=nbreads$NbResampled[which(nbreads$Sample==sample)]
      
      this.reads=read.table(paste(pathExpression, sp, "/", sample,"/ReadCounts_",type,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

      print(length(which(this.reads[,2]>0)))
      
      unique.reads[[sample]]=this.reads[,"NbUniqueReads"]
      names(unique.reads[[sample]])=this.reads[,"GeneID"]
    }
    
    ## reorder values
    
    gene.order=names(unique.reads[[samples[1]]])
    
    for(sample in samples){
      unique.reads[[sample]]=unique.reads[[sample]][gene.order]
    }
    
    ## make data frames
    
    unique.reads=as.data.frame(unique.reads)
    
    ## add gene id as a column
    
    unique.reads$GeneID=gene.order
    unique.reads=unique.reads[,c("GeneID",samples)]
    
    ## write output 
    
    write.table(unique.reads, file=paste(pathExpression, sp, "/AllSamples_UniqueReadCounts_", type,".txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")
    
  }
}

########################################################################
