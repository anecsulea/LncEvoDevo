####################################################################

## get command line arguments

args = commandArgs(trailingOnly=TRUE)

if(length(args)!=5){
  stop("You need to provide 5 arguments: cluster, species, sample, nthreads, strand")
}

cluster=args[1]
species=args[2]
sample=args[3]
nthreads=as.numeric(args[4])
strand=as.numeric(args[5])

####################################################################

if(cluster=="pbil"){
 path="LncEvoDevo/"
} else{
  if(cluster=="in2p3"){
    path="LncEvoDevo/"
  } else{
    stop(paste("Unknown cluster: ",cluster))
  }
}

if(!(species%in%c("Mouse", "Rat", "Chicken"))){
  stop(paste("Unknown species: ",species))
}

pathHisat=paste(path, "results/hisat/", species,"/",sample,"/", sep="")
pathAln=paste(pathHisat, "accepted_hits.bam", sep="")
pathResults=paste(path, "results/expression_estimation/", species, "/",sample, "/", sep="")
pathGeneOverlaps=paste(path, "results/gene_overlaps/", species,"/", sep="")

if(!dir.exists(pathHisat)){
  stop("Weird! cannot find Hisat directory")
}

if(!file.exists(paste(pathAln))){
  stop("Weird! cannot find alignment file")
}

if(!dir.exists(pathResults)){
  stop("Weird! cannot find output directory")
}

release=94

pathAnnot=paste(pathGeneOverlaps, "ExonBlocks_ExcludingOverlapOtherGenes_FilteredTranscripts_StringTie_Ensembl",release,".gtf", sep="")

annot="NonOverlappingExonBlocks_StringTie"

if(!file.exists(pathAnnot)){
  stop(paste("Cannot find annotation file: ",pathAnnot, sep=""))
}

####################################################################

pathOutput=paste(pathResults, "ReadCounts_", annot, ".txt", sep="")

if(file.exists(pathOutput)){
  stop("already done")
}

####################################################################

library(Rsubread)

####################################################################

res.unique=featureCounts(files=pathAln, annot.ext=pathAnnot, isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_id", countMultiMappingReads=FALSE, strandSpecific=strand, nthreads=nthreads)

####################################################################

all.counts=data.frame("GeneID"=rownames(res.unique$counts), "NbUniqueReads"=res.unique$counts[,1], stringsAsFactors=F)

write.table(all.counts, file=paste(pathResults, "ReadCounts_", annot, ".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

####################################################################



