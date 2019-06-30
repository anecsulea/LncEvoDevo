##############################################################################

path="LncEvoDevo/"
pathStringTie=paste(path, "results/stringtie_assembly/", sep="")
pathCufflinks=paste(path, "results/cufflinks_assembly/", sep="")

options(scipen=999) ## remove scientific notation ## options(scipen=0) to get it back

suffixes=list()
suffixes[["Cufflinks"]]=c("merged", "FilteredTranscripts_Cufflinks_Ensembl94")
suffixes[["StringTie"]]=c("assembled_transcripts", "FilteredTranscripts_StringTie_Ensembl94")

paths=list()
paths[["Cufflinks"]]=pathCufflinks
paths[["StringTie"]]=pathStringTie

##############################################################################

for(sp in c("Mouse", "Rat", "Chicken")){

  for(type in c("StringTie")){ ## StringTie already done
    
    pathAssembly=paths[[type]]
      
    for(suffix in suffixes[[type]]){

      print(paste(sp, type, suffix))
      
      cuff=readLines(paste(pathAssembly,sp,"/combined/",suffix,".gtf",sep=""))
      
      print(paste(length(cuff),"exons initially"))
      
      info=lapply(cuff, function(x) unlist(strsplit(x,split="\t")))
      chromo=unlist(lapply(info, function(x) x[1]))
      strand=unlist(lapply(info, function(x) x[7]))
      
      cuff=cuff[which(chromo%in%c(as.character(1:50), "X", "Y", "Z", "W") & strand%in%c("+", "-"))]
      chromo=chromo[which(chromo%in%c(as.character(1:50), "X", "Y", "Z", "W") & strand%in%c("+", "-"))]
      
      print(paste(length(cuff),"exons after selecting chromosomes and strands"))
      
      cuff=paste("chr", cuff,sep="")
      
      header=paste("track type=gtf name='",sp,"'",sep="")
      
      outcuff=c(header, cuff)
      
      writeLines(outcuff, con=paste(pathAssembly,sp,"/combined/",suffix,"_UCSC.gtf",sep=""))
      
      
    }
  }
}

##############################################################################
