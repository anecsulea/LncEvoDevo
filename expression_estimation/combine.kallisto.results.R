########################################################################

path="LncEvoDevo/"
pathExpression=paste(path,"results/expression_estimation/", sep="")
pathEnsembl=paste(path,"data/ensembl_annotations/", sep="")
pathStringTie=paste(path,"results/stringtie_assembly/",sep="")

set.seed(19)

source("normalization.R")

types=c("StringTie")
splist=c("Mouse", "Rat", "Chicken")

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

    read.counts=list()
    tpm=list()
    
    for(sample in samples){

      print(sample)

      this.data=c()
      
      try(this.data<-read.table(paste(pathExpression, sp, "/", sample,"/kallisto_",type, "/abundance.tsv",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\""))

      if(length(this.data)>0){
        this.data$GeneID=unlist(lapply(this.data$target_id, function(x) unlist(strsplit(x, split=":"))[1]))
      
        gene.counts=as.numeric(tapply(as.numeric(this.data$est_counts), as.factor(this.data$GeneID), sum, na.rm=T))
        gene.tpm=as.numeric(tapply(as.numeric(this.data$tpm), as.factor(this.data$GeneID), sum, na.rm=T))

        names(gene.counts)=levels(as.factor(this.data$GeneID))
        names(gene.tpm)=levels(as.factor(this.data$GeneID))
        
        read.counts[[sample]]=gene.counts
        tpm[[sample]]=gene.tpm
      }
    }

    samples=names(read.counts) ## to exclude samples for which we didn't find expression values
    samples.mainstrain=setdiff(samples, samples.otherstrain[[sp]])

    print(paste(length(samples), "samples in total, ",length(samples.mainstrain), "samples for the main strain"))
    
    ## reorder values

    gene.order=names(read.counts[[samples[1]]])
    
    for(sample in samples){
      print(paste("reordering", sample))
      read.counts[[sample]]=read.counts[[sample]][gene.order]
      tpm[[sample]]=tpm[[sample]][gene.order]
    }
    
    ## make data frames
    
    read.counts=as.data.frame(read.counts)
    rownames(read.counts)=gene.order
    
    tpm=as.data.frame(tpm)
    rownames(tpm)=gene.order

    norm.data=normalization(tpm)
    tpm.norm=norm.data[["expdata.norm"]]
    rownames(tpm.norm)=gene.order

    hk.genes=norm.data[["hk.genes"]]
    
    ## add gene id as a column
    
    tpm$GeneID=gene.order
    tpm.norm$GeneID=gene.order
    read.counts$GeneID=gene.order

    tpm=tpm[,c("GeneID", setdiff(colnames(tpm), "GeneID"))]
    tpm.norm=tpm.norm[,c("GeneID", setdiff(colnames(tpm), "GeneID"))]
    read.counts=read.counts[,c("GeneID", setdiff(colnames(tpm), "GeneID"))]

    ## write output 
    
    write.table(read.counts, file=paste(pathExpression, sp, "/AllSamples_KallistoEstimatedCounts_", type,".txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")
    
    write.table(tpm, file=paste(pathExpression, sp, "/AllSamples_KallistoTPM_", type,".txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")

    write.table(tpm.norm, file=paste(pathExpression, sp, "/AllSamples_KallistoNormalizedTPM_", type,".txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")

    ## write output only for main strain samples
    
    write.table(read.counts[,c("GeneID", samples.mainstrain)], file=paste(pathExpression, sp, "/AllSamples_KallistoEstimatedCounts_", type,"_MainStrain.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")

    write.table(tpm[,c("GeneID", samples.mainstrain)], file=paste(pathExpression, sp, "/AllSamples_KallistoTPM_", type,"_MainStrain.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")

    write.table(tpm.norm[,c("GeneID", samples.mainstrain)], file=paste(pathExpression, sp, "/AllSamples_KallistoNormalizedTPM_", type,"_MainStrain.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")
     
  }
}

########################################################################
