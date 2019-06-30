##########################################################################

path="LncEvoDevo/"
pathInfo=paste(path, "results/lncRNA_dataset/", sep="")
pathExpression=paste(path, "results/expression_estimation/", sep="")
pathAnnot=paste(path, "data/ensembl_annotations/", sep="")
pathResults=paste(path, "supplementary_datasets/SupplementaryDataset1/", sep="")

options(scipen=999) ## remove scientific notation ## options(scipen=0) to get it back

release=94

minreads=10 ## we also select only protein-coding genes with a minimum number of reads

##########################################################################

for(sp in c("Mouse", "Rat", "Chicken")){
  genenames=read.table(paste(pathAnnot, sp, "/GeneNames_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(genenames)=genenames[,1]
  genenames=genenames[which(genenames[,2]!=""),]

  allinfo=read.table(paste(pathInfo, sp, "/AllInfo_WithRejectionReason_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

  class.cols=c("Class", grep("BlastX", colnames(allinfo), value=T), grep("CSF", colnames(allinfo), value=T))

  allinfo$GeneName=rep("NA", dim(allinfo)[1])
  allinfo$GeneName[which(allinfo$GeneID%in%rownames(genenames))]=genenames[allinfo$GeneID[which(allinfo$GeneID%in%rownames(genenames))],2]

  if(sp=="Mouse"){
    allinfo=allinfo[, c("GeneID", "GeneName", "Chr", "Start", "End", "Strand", "ExonicLength", "NbExons", "AnnotationSource",  "EnsemblBiotype", class.cols, "UnmappableLength", "AllSupportedJunctions", "ReadCountMainStrain","MaxTPM", "MaxRatioSenseAntisense", "Distance5PC", "Distance3PC", "RegionPCSense", "RegionPCAntisense", "ExpressionCorrelationNoOverlap",  "OverlapRNARepeats", "OverlapRetrogenes", "OverlapTSSEncodeEnhancer", "OverlapTSSVistaEnhancer", "tRNAPrecursor", "BidirectionalPromoter1kb", "BidirectionalPromoterProteinCoding1kb", "CGI500bp", "RejectionReason")]

    colnames(allinfo)[which(colnames(allinfo)=="OverlapTSSEncodeEnhancer")]="OverlapEncodeEnhancer"
    colnames(allinfo)[which(colnames(allinfo)=="OverlapTSSVistaEnhancer")]="OverlapVistaEnhancer"
    
  } else{
    allinfo=allinfo[, c("GeneID", "GeneName", "Chr", "Start", "End", "Strand", "ExonicLength", "NbExons", "AnnotationSource",  "EnsemblBiotype", class.cols, "UnmappableLength", "AllSupportedJunctions", "ReadCountMainStrain","MaxTPM", "MaxRatioSenseAntisense", "Distance5PC", "Distance3PC", "RegionPCSense", "RegionPCAntisense",  "ExpressionCorrelationNoOverlap", "OverlapRNARepeats", "OverlapRetrogenes", "tRNAPrecursor", "BidirectionalPromoter1kb", "BidirectionalPromoterProteinCoding1kb", "CGI500bp","RejectionReason")]
  }

  ## counts for all samples also

  expression=read.table(paste(pathExpression, sp, "/AllSamples_UniqueReadCounts_StringTie.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  rownames(expression)=expression$GeneID
  expression=expression[,-which(colnames(expression)%in%c("GeneID", "ExonicLength"))]

  totalcount=apply(expression, 1, sum)
  names(totalcount)=rownames(expression)

  allinfo$ReadCountAllSamples=totalcount[allinfo$GeneID]
  
  colnames(allinfo)[which(colnames(allinfo)=="tRNAPrecursor")]="OverlaptRNA"
  colnames(allinfo)[which(colnames(allinfo)=="AllSupportedJunctions")]="VerifiedSpliceJunctions"
  colnames(allinfo)[which(colnames(allinfo)=="Distance5PC")]="Distance5ProteinCoding"
  colnames(allinfo)[which(colnames(allinfo)=="Distance3PC")]="Distance3ProteinCoding"
  colnames(allinfo)[which(colnames(allinfo)=="ReadCountMainStrain")]="ReadCount"

  colnames(allinfo)[which(colnames(allinfo)=="BidirectionalPromoter1kb")]="BidirectionalPromoter"
  colnames(allinfo)[which(colnames(allinfo)=="BidirectionalPromoterProteinCoding1kb")]="BidirectionalPromoterProteinCoding"
  colnames(allinfo)[which(colnames(allinfo)=="CGI500bp")]="CGIPromoter"

  colnames(allinfo)[which(colnames(allinfo)=="ExpressionCorrelationNoOverlap")]="ExpCorrNonOverlappingExons"
  
  allinfo=allinfo[,c(setdiff(colnames(allinfo), "RejectionReason"), "RejectionReason")]

  allinfo$CandidateLncRNA=rep("No", dim(allinfo)[1])
  allinfo$CandidateLncRNA[which(allinfo$RejectionReason=="")]="Yes"

  
  allinfo$SelectedProteinCoding=rep("No", dim(allinfo)[1])
  allinfo$SelectedProteinCoding[which(allinfo$EnsemblBiotype=="protein_coding" & allinfo$ReadCount>=minreads)]="Yes"

   
  write.table(allinfo, file=paste(pathResults, "/GeneInfo_",sp, ".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

}

##########################################################################


