########################################################################

path="LncEvoDevo/"
pathExpression=paste(path,"results/expression_estimation/", sep="")
pathEnsembl=paste(path,"data/ensembl_annotations/", sep="")
pathStringTie=paste(path,"results/stringtie_assembly/",sep="")

set.seed(19)

source("normalization.R")

types=c("Ensembl", "StringTie")
#splist=c("Mouse", "Rat", "Chicken")
splist=c("Chicken")

release=94

########################################################################

samples.otherstrain=list()
samples.otherstrain[["Mouse"]]=c("Sertoli", "Spermatids", "Spermatocytes", "Spermatogonia", "Spermatozoa")
samples.otherstrain[["Rat"]]=c("Brain_Adult3", "Brain_Adult4", "Kidney_Adult3", "Kidney_Adult4", "Liver_Adult3", "Liver_Adult4", "Testis_Adult4")
samples.otherstrain[["Chicken"]]=c("Brain_Adult1", "Brain_Adult2", "Brain_Day18Embryo1", "Brain_Day18Embryo10", "Brain_Day18Embryo2", "Brain_Day18Embryo3", "Brain_Day18Embryo4", "Brain_Day18Embryo5", "Brain_Day18Embryo6", "Brain_Day18Embryo7", "Brain_Day18Embryo8", "Brain_Day18Embryo9",  "Kidney_Adult1", "Kidney_Adult2", "Kidney_Day18Embryo1", "Kidney_Day18Embryo10", "Kidney_Day18Embryo2", "Kidney_Day18Embryo3", "Kidney_Day18Embryo4", "Kidney_Day18Embryo5", "Kidney_Day18Embryo6", "Kidney_Day18Embryo7", "Kidney_Day18Embryo8", "Kidney_Day18Embryo9", "Liver_Adult1", "Liver_Adult2", "Liver_Day18Embryo1", "Liver_Day18Embryo3", "Liver_Day18Embryo4", "Liver_Day18Embryo5", "Liver_Day18Embryo6", "Liver_Day18Embryo7", "Liver_Day18Embryo8", "Liver_Day18Embryo9", "Testis_Adult1", "Testis_Day18Embryo1", "Testis_Day18Embryo2", "Testis_Day18Embryo3", "Testis_Day18Embryo4", "Testis_Day18Embryo5", "Testis_Day4.5Embryo1", "Testis_Day4.5Embryo2", "Testis_Day6Embryo1", "Testis_Day6Embryo2")

########################################################################

for(sp in splist){
  print(sp)
  
  samples=system(paste("ls ", pathExpression, sp, " | grep -v txt", sep=""), intern=T)

  for(type in types){
    print(type)

    if(type=="Ensembl"){
      exons=read.table(paste(pathEnsembl,sp,"/ExonBlocks_FilteredTranscripts_Ensembl",release,".txt",sep=""), h=F, stringsAsFactors=F, sep="\t")
    }

    if(type=="StringTie"){
       exons=read.table(paste(pathStringTie,sp,"/combined/ExonBlocks_FilteredTranscripts_StringTie_Ensembl",release,".txt",sep=""), h=F, stringsAsFactors=F, sep="\t")
    }
    
    exons$Length=exons$V5-exons$V4+1
    
    exonic.length=as.numeric(tapply(exons$Length, as.factor(exons$V1), sum))
    names(exonic.length)=levels(as.factor(exons$V1))
    
    
    unique.reads=list()
     
    for(sample in samples){

      print(sample)

      this.reads=c()
      
      try(this.reads<-read.table(paste(pathExpression, sp, "/", sample,"/ReadCounts_",type,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\""))

      if(length(this.reads)>0){
      ## remove duplicates
        unique.reads[[sample]]=this.reads[,"NbUniqueReads"]
        names(unique.reads[[sample]])=this.reads[,"GeneID"]
      }
    }

    samples=names(unique.reads) ## to exclude samples for which we didn't find expression values
    samples.mainstrain=setdiff(samples, samples.otherstrain[[sp]])

    print(paste(length(samples), "samples in total, ",length(samples.mainstrain), "samples for the main strain"))
    
    ## reorder values

    exonic.length=exonic.length[intersect(names(exonic.length), names(unique.reads[[samples[1]]]))]
    
    gene.order=names(exonic.length)
    
    for(sample in samples){
      print(paste("reordering", sample))
      unique.reads[[sample]]=unique.reads[[sample]][gene.order]
    }
    
    ## make data frames
    
    unique.reads=as.data.frame(unique.reads)
    
    ## add gene id as a column
    
    unique.reads$GeneID=gene.order
    unique.reads$ExonicLength=exonic.length
    unique.reads=unique.reads[,c("GeneID","ExonicLength", samples)]

    ## write output 
    
    write.table(unique.reads, file=paste(pathExpression, sp, "/AllSamples_UniqueReadCounts_", type,".txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")

    ## write output only for main strain samples
    
    unique.reads.mainstrain=unique.reads[,c("GeneID","ExonicLength", samples.mainstrain)]
    
    write.table(unique.reads.mainstrain, file=paste(pathExpression, sp, "/AllSamples_UniqueReadCounts_", type,"_MainStrain.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")
     
  }
}

########################################################################
