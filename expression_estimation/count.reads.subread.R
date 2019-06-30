####################################################################

## get command line arguments

args = commandArgs(trailingOnly=TRUE)

if(length(args)!=7){
  stop("You need to provide 7 arguments: cluster, annotation, Ensembl release number, species, sample, nthreads, strand")
}

cluster=args[1]
annot=args[2]
release=args[3]
species=args[4]
sample=args[5]
nthreads=as.numeric(args[6])
strand=as.numeric(args[7])

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
pathEnsembl=paste(path, "data/ensembl_annotations/", species,"/", sep="")
pathStringTie=paste(path, "results/stringtie_assembly/", species, "/combined/", sep="")

if(!dir.exists(pathHisat)){
  stop("Weird! cannot find Hisat directory")
}

if(!file.exists(paste(pathAln))){
  stop("Weird! cannot find alignment file")
}

if(!dir.exists(pathResults)){
  stop("Weird! cannot find output directory")
}

if(annot=="Ensembl"){
  pathAnnot=paste(pathEnsembl, "FilteredTranscripts_Ensembl",release,".gtf", sep="")
} else{
  if(annot=="StringTie"){
    pathAnnot=paste(pathStringTie, "FilteredTranscripts_StringTie_Ensembl", release, ".gtf", sep="")
  }
}

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

res.all=featureCounts(files=pathAln, annot.ext=pathAnnot, isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_id", countMultiMappingReads=TRUE, strandSpecific=strand, nthreads=nthreads)

if(!(all(rownames(res.unique$counts)%in%rownames(res.all$counts)))){
  stop("Weird! not the same genes in unique and non-unique counts")
}

if(!(all(rownames(res.all$counts)%in%rownames(res.unique$counts)))){
  stop("Weird! not the same genes in unique and non-unique counts")
}

####################################################################

all.counts=data.frame("GeneID"=rownames(res.unique$counts), "NbUniqueReads"=res.unique$counts[,1], stringsAsFactors=F)

all.counts$NbReadsTotal=res.all$counts[all.counts$GeneID, 1]

write.table(all.counts, file=paste(pathResults, "ReadCounts_", annot, ".txt", sep=""), row.names=F, col.names=T, sep="\t")

####################################################################



