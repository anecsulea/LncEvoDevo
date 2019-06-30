########################################################################

path="LncEvoDevo/"
pathHisat=paste(path,"results/hisat/", sep="")
pathResults=paste(path,"results/antisense_transcription/", sep="")
pathAnnot=paste(path,"data/ensembl_annotations/", sep="")
pathStringTie=paste(path,"results/stringtie_assembly/",sep="")

set.seed(19)


types=c("FilteredTranscripts_StringTie_Ensembl94")

########################################################################

for(sp in c("Mouse", "Rat", "Chicken")){
  print(sp)
  
  samples=system(paste("ls ", pathHisat, sp," | grep -v txt", sep=""), intern=T)
  
  for(type in types){
    print(type)
    
    coverage=list()
    info=c()
     
    for(sample in samples){

      print(sample)
      
      this.cov=read.table(paste(pathResults, sp, "/", sample,"/CoverageTranscripts_",type,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
      
      if(length(info)==0){
        info=this.cov[,c("GeneID", "TranscriptID")]
        rownames(info)=info$TranscriptID
      }
      
      ## remove duplicates
      
      coverage[[sample]]=this.cov[,c("CoverageSense", "CoverageAntisense")]
      rownames(coverage[[sample]])=this.cov[,"TranscriptID"]
    }
    
    ## reorder values
    
    transcript.order=rownames(coverage[[samples[1]]])
    
    for(sample in samples){
      print(paste("reordering", sample))
      coverage[[sample]]=coverage[[sample]][transcript.order,]
    }
    
    ## make data frames
    
    coverage=as.data.frame(coverage)
    
    ## add gene id as a column
    
    coverage$TranscriptID=transcript.order
    coverage$GeneID=info[transcript.order, "GeneID"]
    
    coverage=coverage[,c("GeneID", "TranscriptID", setdiff(colnames(coverage), c("GeneID", "TranscriptID")))]
    

    ## write output 
    
    write.table(coverage, file=paste(pathResults, sp, "/AllSamples_CoverageTranscripts_", type,".txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")
  }
}

########################################################################
