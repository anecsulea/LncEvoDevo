#############################################################################

path="LncEvoDevo/"

pathResults=paste(path, "results/lncRNA_dataset/",sep="")

release=94

mindist.pc=5000 ## for de novo-annotated lncRNAs
minreads=10
maxfrrna=0.25 ## we don't tolerate large overlaps with RNA repeats
maxfrretro=0.5 ## we don't tolerate large (>50%) overlaps with retrogenes
maxfrallrep=1 ## no constraint on repeats in general
minlength.multiex=200
minlength.monoex=500
minratio.sas=0.01
maxfrunmap=0.05
minexpcor=0.9

acceptable.biotypes=c("lincRNA", "processed_transcript", "antisense", "TEC", "bidirectional_promoter_lncRNA", "macro_lncRNA", "sense_intronic")

#############################################################################

for(sp in c("Mouse", "Rat", "Chicken")){

  print(sp)

  allinfo=read.table(paste(pathResults, sp, "/AllInfo_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  rownames(allinfo)=allinfo$GeneID

  allinfo$FractionUnmappable=allinfo$UnmappableLength/allinfo$ExonicLength

  allinfo$RejectionReason=rep("", dim(allinfo)[1])
  
  ## select lncRNAs

  ## acceptable biotypes
  
  all.lnc=allinfo$GeneID[which(is.na(allinfo$EnsemblBiotype) | allinfo$EnsemblBiotype%in%acceptable.biotypes)]
  print(paste(length(all.lnc), "de novo transcripts or Ensembl with acceptable biotypes"))

  notok=which((!is.na(allinfo$EnsemblBiotype)) & (!(allinfo$EnsemblBiotype%in%acceptable.biotypes)))
  allinfo$RejectionReason[notok]=paste("bad_biotype", allinfo$RejectionReason[notok], sep=",")
  
  ## classified as noncoding, for new genes
  
  all.lnc=all.lnc[which((allinfo[all.lnc, "Class"]=="noncoding") | (allinfo[all.lnc, "AnnotationSource"]=="Ensembl"))]
  print(paste(length(all.lnc), "noncoding transcripts"))
  
  notok=which((allinfo[, "Class"]=="coding" | is.na(allinfo[, "Class"])) & allinfo[,"AnnotationSource"]=="DeNovo")
  allinfo$RejectionReason[notok]=paste("coding_potential", allinfo$RejectionReason[notok], sep=",")
  
  ## exonic length, only for new genes
  
  all.lnc=all.lnc[which((allinfo[all.lnc,"ExonicLength"]>=minlength.multiex & allinfo[all.lnc,"NbExons"]>=2) | (allinfo[all.lnc,"ExonicLength"]>=minlength.monoex & allinfo[all.lnc,"NbExons"]==1) | (allinfo[all.lnc, "AnnotationSource"] == "Ensembl")) ]
  print(paste(length(all.lnc), "long transcripts"))

  notok=which(allinfo[,"ExonicLength"]<minlength.monoex & allinfo[,"NbExons"]==1 & allinfo[, "AnnotationSource"]=="DeNovo")
  allinfo$RejectionReason[notok]=paste("too_short_monoexonic", allinfo$RejectionReason[notok], sep=",")

  notok=which(allinfo[,"ExonicLength"]<minlength.multiex & allinfo[,"NbExons"]>=2  & allinfo[, "AnnotationSource"]=="DeNovo")
  allinfo$RejectionReason[notok]=paste("too_short_multiexonic", allinfo$RejectionReason[notok], sep=",")

  ## fraction unmappable, for new genes
  
  all.lnc=all.lnc[which((allinfo[all.lnc,"FractionUnmappable"]<=maxfrunmap) | (allinfo[all.lnc, "AnnotationSource"] == "Ensembl"))]
  print(paste(length(all.lnc), "transcripts with at most", maxfrunmap," unmappable exonic length"))

  notok=which(allinfo[,"FractionUnmappable"]>maxfrunmap & allinfo[,"AnnotationSource"]=="DeNovo")
  allinfo$RejectionReason[notok]=paste("fraction_unmappable", allinfo$RejectionReason[notok], sep=",")

  ## region PC sense  - for de novo genes

  all.lnc=all.lnc[which((allinfo[all.lnc,"RegionPCSense"]=="intergenic") | (allinfo[all.lnc, "AnnotationSource"] == "Ensembl"))]
  print(paste(length(all.lnc), "intergenic"))

  notok=which(allinfo[,"RegionPCSense"]!="intergenic" &  allinfo[, "AnnotationSource"]=="DeNovo")
  allinfo$RejectionReason[notok]=paste("region_PC_sense", allinfo$RejectionReason[notok], sep=",")


  ## far from protein-coding gene exons, for de novo genes
  
  all.lnc=all.lnc[which((allinfo[all.lnc,"Distance5PC"]>=mindist.pc & allinfo[all.lnc,"Distance3PC"]>=mindist.pc) | (allinfo[all.lnc, "AnnotationSource"]=="Ensembl"))]
  print(paste(length(all.lnc), "far from protein-coding gene exons"))

  notok=which((allinfo[,"Distance5PC"]<mindist.pc | allinfo[,"Distance3PC"]<mindist.pc) & allinfo[, "AnnotationSource"]=="DeNovo")
  allinfo$RejectionReason[notok]=paste("distance_PC", allinfo$RejectionReason[notok], sep=",")

  ## ratio sense-antisense transcription, for all genes

  all.lnc=all.lnc[which(allinfo[all.lnc,"MaxRatioSenseAntisense"]>=minratio.sas)]
  print(paste(length(all.lnc), "with at least",minratio.sas,"sense-antisense transcription in at least one sample"))
  
  notok=which(allinfo[,"MaxRatioSenseAntisense"]<minratio.sas | is.na(allinfo[,"MaxRatioSenseAntisense"]))
  allinfo$RejectionReason[notok]=paste("ratio_sense_antisense", allinfo$RejectionReason[notok], sep=",")

  ## min number of reads - across our samples only, because for chicken we don't control that
 
  all.lnc=all.lnc[which(allinfo[all.lnc,"ReadCountMainStrain"]>=minreads)]
  print(paste(length(all.lnc), "supported by at least", minreads, "reads"))

  notok=which(allinfo[,"ReadCountMainStrain"]<minreads)
  allinfo$RejectionReason[notok]=paste("read_count", allinfo$RejectionReason[notok], sep=",")
  
  ## overlap RNA repeats

  all.lnc=all.lnc[which(allinfo[all.lnc,"OverlapRNARepeats"]<=maxfrrna)]
  print(paste(length(all.lnc), "with less than", maxfrrna, "overlap with RNA repeats"))

  notok=which(allinfo[,"OverlapRNARepeats"]>maxfrrna)
  allinfo$RejectionReason[notok]=paste("overlap_RNA_repeats", allinfo$RejectionReason[notok], sep=",")

  ## overlap all repeats
  all.lnc=all.lnc[which(allinfo[all.lnc,"OverlapAllRepeats"]<=maxfrallrep)]
  print(paste(length(all.lnc), "with less than", maxfrallrep, "overlap with all types of repeats"))

  notok=which(allinfo[,"OverlapAllRepeats"]>maxfrallrep)
  allinfo$RejectionReason[notok]=paste("overlap_all_repeats", allinfo$RejectionReason[notok], sep=",")

  ## overlap retrogenes

  all.lnc=all.lnc[which(allinfo[all.lnc,"OverlapRetrogenes"]<=maxfrretro)]
  print(paste(length(all.lnc), "with less than", maxfrretro, "overlap with retrogenes"))
  
  notok=which(allinfo[,"OverlapRetrogenes"]>maxfrretro)
  allinfo$RejectionReason[notok]=paste("overlap_retrogenes", allinfo$RejectionReason[notok], sep=",")
  

  ## tRNA precursors
  all.lnc=all.lnc[which(allinfo[all.lnc,"tRNAPrecursor"]=="No")]
  print(paste(length(all.lnc), "that are not tRNA precursors"))

  notok=which(allinfo[,"tRNAPrecursor"]!="No")
  allinfo$RejectionReason[notok]=paste("tRNA", allinfo$RejectionReason[notok], sep=",")

  ## supported junctions
  
  all.lnc=all.lnc[which((allinfo[all.lnc,"AnnotationSource"]=="DeNovo" & allinfo[all.lnc, "HasWrongStrandJunctions"]=="No") | (allinfo[all.lnc,"AnnotationSource"]=="Ensembl"))]
  print(paste(length(all.lnc), "that are either Ensembl-annotated, or de novo without wrong strand junctions"))

  notok=which(allinfo[,"AnnotationSource"]=="DeNovo" & allinfo[, "HasWrongStrandJunctions"]=="Yes")
  allinfo$RejectionReason[notok]=paste("not_supported_junctions", allinfo$RejectionReason[notok], sep=",")

  ## expression correlation between all regions and no exonic overlap

  all.lnc=all.lnc[which(allinfo[all.lnc,"ExpressionCorrelationNoOverlap"]>=minexpcor)]
  print(paste(length(all.lnc), "that have good expression correlation between all regions and no-overlapping regions"))
  
  notok=which(allinfo[,"ExpressionCorrelationNoOverlap"]<minexpcor | is.na(allinfo[,"ExpressionCorrelationNoOverlap"]))
  allinfo$RejectionReason[notok]=paste("expression_correlation_nooverlap", allinfo$RejectionReason[notok], sep=",")
  

  ## check if everything is ok

  if(length(which(allinfo$RejectionReason==""))!=length(all.lnc)){
    stop("Weird ! discrepancies between rejection reason and lncRNA numbers")
  }
  
  if(length(which(allinfo[all.lnc,"RejectionReason"]==""))!=length(all.lnc)){
    stop("Weird ! some lncRNAs should have been rejected")
  }

  print(sort(table(allinfo$RejectionReason)))
  
  print(paste(length(all.lnc), "lncRNAs in total"))
  results=allinfo[all.lnc,]
  
  print(table(results$RegionPCSense))
  print(table(results$RegionPCAntisense))
  print(table(results$AnnotationSource))
  print(table(results$NbExons))
  print(table(results$EnsemblBiotype))
    
  ## write results

  write.table(results, file=paste(pathResults, sp, "/SelectedLncRNAs_FilteredTranscripts_StringTie_Ensembl", release, ".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

  ## write results + rejection reason

  write.table(allinfo, file=paste(pathResults, sp, "/AllInfo_WithRejectionReason_FilteredTranscripts_StringTie_Ensembl", release, ".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
}

#############################################################################

