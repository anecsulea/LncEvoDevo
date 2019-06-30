#############################################################################

path="LncEvoDevo/"
pathStringTie=paste(path, "results/stringtie_assembly/", sep="")
pathAntisenseTranscription=paste(path, "results/antisense_transcription/", sep="")
pathCodingPotential=paste(path, "results/coding_potential/", sep="")
pathGeneOverlaps=paste(path, "results/gene_overlaps/",sep="")
pathEnsembl=paste(path, "data/ensembl_annotations/", sep="")
pathExpression=paste(path, "results/expression_estimation/", sep="")
pathRepeats=paste(path, "results/overlap_repeats/", sep="")
pathUnmap=paste(path, "results/overlap_unmappable_regions/", sep="")
pathOrtho=paste(path, "results/ortho_genes/whole_genome_alignments/", sep="")
pathPCClusters=paste(path, "results/protein_coding_clusters/", sep="")
pathEnhancers=paste(path, "results/overlap_enhancers/", sep="")
pathResults=paste(path, "results/lncRNA_dataset/",sep="")

#############################################################################

release=94

#############################################################################

mindist.pc=5000
minreads=10
minratio.sas=0.01
minfrrna=0.1

splist=c("Mouse", "Rat", "Chicken")
samples.otherstrain=list()
samples.otherstrain[["Mouse"]]=c("Sertoli", "Spermatids", "Spermatocytes", "Spermatogonia", "Spermatozoa")
samples.otherstrain[["Rat"]]=c("Brain_Adult3", "Brain_Adult4", "Kidney_Adult3", "Kidney_Adult4", "Liver_Adult3", "Liver_Adult4", "Testis_Adult4")
samples.otherstrain[["Chicken"]]=c("Brain_Adult1", "Brain_Adult2", "Brain_Day18Embryo1", "Brain_Day18Embryo10", "Brain_Day18Embryo2", "Brain_Day18Embryo3", "Brain_Day18Embryo4", "Brain_Day18Embryo5", "Brain_Day18Embryo6", "Brain_Day18Embryo7", "Brain_Day18Embryo8", "Brain_Day18Embryo9",  "Kidney_Adult1", "Kidney_Adult2", "Kidney_Day18Embryo1", "Kidney_Day18Embryo10", "Kidney_Day18Embryo2", "Kidney_Day18Embryo3", "Kidney_Day18Embryo4", "Kidney_Day18Embryo5", "Kidney_Day18Embryo6", "Kidney_Day18Embryo7", "Kidney_Day18Embryo8", "Kidney_Day18Embryo9", "Liver_Adult1", "Liver_Adult2", "Liver_Day18Embryo1", "Liver_Day18Embryo3", "Liver_Day18Embryo4", "Liver_Day18Embryo5", "Liver_Day18Embryo6", "Liver_Day18Embryo7", "Liver_Day18Embryo8", "Liver_Day18Embryo9", "Testis_Adult1", "Testis_Day18Embryo1", "Testis_Day18Embryo2", "Testis_Day18Embryo3", "Testis_Day18Embryo4", "Testis_Day18Embryo5", "Testis_Day4.5Embryo1", "Testis_Day4.5Embryo2", "Testis_Day6Embryo1", "Testis_Day6Embryo2")

#############################################################################

for(sp in c("Mouse", "Rat", "Chicken")){

  print(sp)

  all.samples=system(paste("ls ", pathStringTie, sp, "/ | grep -v combined", sep=""), intern=T)
  samples=setdiff(all.samples, samples.otherstrain[[sp]])

  print(paste(length(samples), "samples for", sp))

  ## read gene info from Ensembl
  
  annot=read.table(paste(pathEnsembl, sp, "/GeneInfo_Ensembl", release, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(annot)=annot$stable_id
  ensembl.genes=annot$stable_id

  ens.nc=annot$stable_id[which(annot$biotype%in%c("lincRNA", "processed_transcript", "antisense", "bidirectional_promoter_lncRNA", "macro_lncRNA", "3prime_overlapping_ncRNA"))]
  ens.unclear=annot$stable_id[which(annot$biotype%in%c("TEC"))]
  ens.pc=annot$stable_id[which(annot$biotype%in%c("protein_coding"))]


  ## read exon blocks for exonic length

  blocks=read.table(paste(pathStringTie, sp, "/combined/ExonBlocks_FilteredTranscripts_StringTie_Ensembl", release, ".txt", sep=""), h=F, stringsAsFactors=F, sep="\t")
  colnames(blocks)=c("GeneID", "ExonID", "Chr", "Start", "End", "Strand")
  
  blocks$Length=blocks$End-blocks$Start+1

  exoniclength=tapply(blocks$Length, as.factor(blocks$GeneID), sum)
  names(exoniclength)=levels(as.factor(blocks$GeneID))

  nbexons=as.numeric(table(blocks$GeneID))
  names(nbexons)=levels(as.factor(blocks$GeneID))

  ## read coding potential

  class=read.table(paste(pathCodingPotential, sp, "/CombinedGeneClassification_FilteredTranscripts_StringTie_Ensembl", release, ".txt", sep=""),  h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(class)=class$GeneID

  ## distance to pc genes

  distpc.ens=read.table(paste(pathGeneOverlaps, sp, "/Distance_ProteinCodingGeneExons_SenseStrand_AllTranscripts_Ensembl", release, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  rownames(distpc.ens)=distpc.ens$GeneID
  far.ens=distpc.ens$GeneID[which(distpc.ens$DistanceLeft>=mindist.pc & distpc.ens$DistanceRight>=mindist.pc)]
  

  distpc.st=read.table(paste(pathGeneOverlaps, sp, "/Distance_ProteinCodingGeneExons_SenseStrand_FilteredTranscripts_StringTie_Ensembl", release, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  rownames(distpc.st)=distpc.st$GeneID
  far.st=distpc.st$GeneID[which(distpc.st$DistanceLeft>=mindist.pc & distpc.st$DistanceRight>=mindist.pc)]

  distpc.st=distpc.st[rownames(distpc.ens),]
  alldist=data.frame("GeneID"=rownames(distpc.ens), "Chr"=distpc.ens$Chr, "Start"=distpc.ens$Start, "End"=distpc.ens$End, "Strand"=distpc.ens$Strand, "DistanceLeft.Ensembl"=distpc.ens$DistanceLeft, "DistanceRight.Ensembl"=distpc.ens$DistanceRight, "DistanceLeft.StringTie"=distpc.st$DistanceLeft, "DistanceRight.StringTie"=distpc.st$DistanceRight)
  alldist$MinDistanceLeft=apply(alldist[,c("DistanceLeft.Ensembl", "DistanceLeft.StringTie")], 1, min)
  alldist$MinDistanceRight=apply(alldist[,c("DistanceRight.Ensembl", "DistanceRight.StringTie")], 1, min)

  rownames(alldist)=alldist$GeneID

  alldist$Distance5=alldist$MinDistanceLeft
  alldist$Distance5[which(alldist$Strand==-1)]=alldist$MinDistanceRight[which(alldist$Strand==-1)]
  
  alldist$Distance3=alldist$MinDistanceRight
  alldist$Distance3[which(alldist$Strand==-1)]=alldist$MinDistanceLeft[which(alldist$Strand==-1)]
  
  ## overlap RNA repeats

  ovrnarep=read.table(paste(pathRepeats, sp, "/OverlapRepeats_RNA_BothStrands_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F)

  rownames(ovrnarep)=ovrnarep$GeneID
  ovrnarep$TotalFraction=apply(ovrnarep[,which(colnames(ovrnarep)!="GeneID")],1, sum)

  ## overlap all repeats

  ovallrep=read.table(paste(pathRepeats, sp, "/OverlapRepeats_TEFamily_BothStrands_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F)
  rownames(ovallrep)=ovallrep$GeneID
  ovallrep$TotalFraction=apply(ovallrep[,which(colnames(ovallrep)!="GeneID")],1, sum)

  ## overlap retrogenes

  ovretro=read.table(paste(pathRepeats, sp, "/OverlapRetrogenes_BothStrands_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F)
  rownames(ovretro)=ovretro$GeneID
  ovretro$TotalFraction=ovretro$OverlapLength/ovretro$TotalLength

  
  ## region with respect to protein-coding genes, on the same strand - we do not keep intronic transcripts

  regions=read.table(paste(pathGeneOverlaps, sp, "/GeneClassification_ProteinCodingGenes_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(regions)=regions$GeneID

  ## sense-antisense transcription
  
  tx.sas=read.table(paste(pathAntisenseTranscription, sp, "/AllSamples_CoverageTranscripts_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, sep="\t", stringsAsFactors=F)
  info.sas=tx.sas[,c("GeneID", "TranscriptID")]
  
  cov.sense=tx.sas[,paste(samples, "CoverageSense", sep=".")]
  cov.antisense=tx.sas[,paste(samples, "CoverageAntisense", sep=".")]
  cov.total=cov.sense+cov.antisense

  ratio=as.matrix(cov.sense)/as.matrix(cov.total)
  maxratio=apply(ratio,1, max, na.rm=T)
  names(maxratio)=tx.sas$TranscriptID
  
   
  maxratio.gene=tapply(maxratio, as.factor(tx.sas$GeneID), max)
  names(maxratio.gene)=levels(as.factor(tx.sas$GeneID))
  
  ## total number of reads

  nbreads=read.table(paste(pathExpression, sp, "/AllSamples_UniqueReadCounts_StringTie.txt", sep=""), h=T, sep="\t", stringsAsFactors=F)
  rownames(nbreads)=nbreads$GeneID

  nbreads=nbreads[,-c(1,2)] ## remove gene id and exonic length info

  print(paste(dim(nbreads)[2], "samples in total"))
  print(paste(length(which(!(colnames(nbreads)%in%samples.otherstrain[[sp]]))), "samples for the main strain"))

  totreads=apply(nbreads, 1, sum)
  names(totreads)=rownames(nbreads)

  totreads.correctstrain=apply(nbreads[,which(!(colnames(nbreads)%in%samples.otherstrain[[sp]]))], 1, sum)
  names(totreads.correctstrain)=rownames(nbreads)

  ## tpm

  tpm=read.table(paste(pathExpression, sp, "/AllSamples_KallistoTPM_StringTie_MainStrain.txt", sep=""),  h=T, sep="\t", stringsAsFactors=F)
  rownames(tpm)=tpm$GeneID

  tpm=tpm[,setdiff(colnames(tpm), "GeneID")]
  
  maxtpm=apply(tpm,1, max)
  names(maxtpm)=rownames(tpm)
  

  ## expression estimation with non-overlapping exonic regions

  expcomp=read.table(paste(pathExpression, sp, "/Comparison_AllRegions_NonOverlappingExonBlocks_StringTie_MainStrain.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  rownames(expcomp)=expcomp$GeneID
  

  ## unmappable regions

  ovunmap=read.table(paste(pathUnmap, sp, "/OverlapUnmappableRegions_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, sep="\t", stringsAsFactors=F)
  rownames(ovunmap)=ovunmap$GeneID

  ## ortho in other species

  for(othersp in setdiff(splist, sp)){
    i=which(splist==sp)
    j=which(splist==othersp)

    if(i<j){
      ortho=read.table(paste(pathOrtho, "ReciprocalBestHits_",sp,"_",othersp,"_StringTie.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
    } else{
      ortho=read.table(paste(pathOrtho, "ReciprocalBestHits_",othersp,"_",sp,"_StringTie.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
    }

    rownames(ortho)=ortho[,paste("ID.",sp,sep="")]

    assign(paste("ortho.",othersp,sep=""),ortho)
  }

  ## antisense TSS

  anti5kb=read.table(paste(pathGeneOverlaps, sp, "/DistanceAntisenseTSS_MaxDist5kb_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F)
  closetss5kb=anti5kb[which(anti5kb$MinDistanceLeft!=">5000" | anti5kb$MaxDistanceRight!=">5000"),]
  closegenes5kb=tapply(closetss5kb$GenesCloseTSS, as.factor(closetss5kb$GeneID), function(x) paste(sort(unique(x)), collapse=","))
  closegenes5kb=unlist(lapply(closegenes5kb, function(x) paste(unique(unlist(strsplit(x,split=","))), collapse=",")))
  names(closegenes5kb)=levels(as.factor(closetss5kb$GeneID))

  closepcgenes5kb=unlist(lapply(closegenes5kb, function(x) paste(intersect(unlist(strsplit(x,split=",")), ens.pc), collapse=",")))
  closepcgenes5kb[which(closepcgenes5kb=="")]=NA
  names(closepcgenes5kb)=names(closegenes5kb)



  anti1kb=read.table(paste(pathGeneOverlaps, sp, "/DistanceAntisenseTSS_MaxDist1kb_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F)
  closetss1kb=anti1kb[which(anti1kb$MinDistanceLeft!=">1000" | anti1kb$MaxDistanceRight!=">1000"),]
  closegenes1kb=tapply(closetss1kb$GenesCloseTSS, as.factor(closetss1kb$GeneID), function(x) paste(sort(unique(x)), collapse=","))
  closegenes1kb=unlist(lapply(closegenes1kb, function(x) paste(unique(unlist(strsplit(x,split=","))), collapse=",")))
  names(closegenes1kb)=levels(as.factor(closetss1kb$GeneID))

  closepcgenes1kb=unlist(lapply(closegenes1kb, function(x) paste(intersect(unlist(strsplit(x,split=",")), ens.pc), collapse=",")))
  closepcgenes1kb[which(closepcgenes1kb=="")]=NA
  names(closepcgenes1kb)=names(closegenes1kb)


  ## overlap with CGI

  cgi=read.table(paste(pathGeneOverlaps, sp, "/OverlapCGI_MaxDist500bp_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F)

  withcgi=cgi[which(!is.na(cgi$OverlappingCGI)),]
  cgigenes=tapply(withcgi$OverlappingCGI, as.factor(withcgi$GeneID), function(x) paste(sort(unique(x)), collapse=";"))
  cgigenes=unlist(lapply(cgigenes, function(x) paste(unique(unlist(strsplit(x,split=";"))), collapse=";")))
  names(cgigenes)=levels(as.factor(withcgi$GeneID))

  ## projected protein-coding clusters

  pc.projection.clusters=read.table(paste(pathPCClusters, sp, "/ProteinCodingGenes_ProjectionClusters_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=F, stringsAsFactors=F, sep="\t")
  colnames(pc.projection.clusters)=c("GeneID","Source")
  
  pc.trinity.clusters=read.table(paste(pathPCClusters, sp, "/ProteinCodingGenes_TrinityClusters_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=F, stringsAsFactors=F, sep="\t")
  colnames(pc.trinity.clusters)=c("GeneID","Source")
  
  
  ## overlap miRNAs

  overlap.mirnas=read.table(paste(pathGeneOverlaps, sp, "/ExonOverlap_StringTie_Ensembl_miRNA.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(overlap.mirnas)=overlap.mirnas$GeneID

  mirna.precursors=overlap.mirnas$GeneID[which(overlap.mirnas$LengthOverlapSense>0)]

  ## overlap tRNAs

  trna.precursors=c()
  
  if(sp=="Mouse"){
    overlap.trnas=read.table(paste(pathGeneOverlaps, sp, "/Overlap_tRNAs_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    rownames(overlap.trnas)=overlap.trnas$GeneID

    trna.precursors=overlap.trnas$GeneID[which(overlap.trnas$LengthOverlapSense>0)]
  }

  ## splice junctions statistics

  stats.junc=read.table(paste(pathStringTie, sp, "/combined/FilteredTranscripts_StringTie_Ensembl",release,"_SpliceJunctionsStats.txt", sep=""), h=T, stringsAsFactors=F)

  stats.junc$FractionSupported=rep(1, dim(stats.junc)[1])
  stats.junc$FractionSupported[which(stats.junc$NbIntrons>0)]=stats.junc$NbSupportedIntrons[which(stats.junc$NbIntrons>0)]/stats.junc$NbIntrons[which(stats.junc$NbIntrons>0)]
  
  nb.tx.wrong.junctions=tapply(stats.junc$NbWrongStrand, as.factor(stats.junc$GeneID), function(x) length(which(x>0)))
  nb.tx.supported=tapply(stats.junc$FractionSupported, as.factor(stats.junc$GeneID), function(x) length(which(x==1)))
  nb.tx.all=table(as.factor(stats.junc$GeneID))

  if(!(all(names(nb.tx.all)==names(nb.tx.supported)))){
    stop("weird transcript names!")
  }
  
  genes.wrong.junctions=names(nb.tx.wrong.junctions)[which(nb.tx.wrong.junctions>0)]
  genes.all.supported=names(nb.tx.all)[which(nb.tx.all==nb.tx.supported)]
  genes.some.supported=names(nb.tx.supported)[which(nb.tx.supported>0)]

  
  all.genes=unique(blocks$GeneID)

  ## overlap between TSS and enhancers

  if(sp=="Mouse"){
    data.oven1000.encode=read.table(paste(pathEnhancers, sp, "/OverlapTSS_EncodeYueLabEnhancers_MaxDist1kb_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    oven1000.encode=rep("No", length(all.genes))
    names(oven1000.encode)=all.genes
    oven1000.encode[which(all.genes%in%data.oven1000.encode$GeneID[which(!is.na(data.oven1000.encode$OverlappingEnhancers))])]="Yes"

    data.oven1000.vista=read.table(paste(pathEnhancers, sp, "/OverlapTSS_Vista_MaxDist1kb_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    oven1000.Vista=rep("No", length(all.genes))
    names(oven1000.Vista)=all.genes
    oven1000.Vista[which(all.genes%in%data.oven1000.vista$GeneID[which(!is.na(data.oven1000.vista$OverlappingEnhancers))])]="Yes"

  } else{
    oven1000.encode=rep(NA, length(all.genes))
    names(oven1000.encode)=all.genes
    
    oven1000.Vista=rep(NA, length(all.genes))
    names(oven1000.Vista)=all.genes
  }
  
    
  ## combine all data

  results=data.frame("GeneID"=all.genes, stringsAsFactors=F)

  results$AnnotationSource=rep("DeNovo", dim(results)[1])
  results$AnnotationSource[which(results$GeneID%in%ensembl.genes)]="Ensembl"
  results$EnsemblBiotype=annot[results$GeneID,"biotype"]
  
  results$Chr=alldist[all.genes,"Chr"]
  results$Start=alldist[all.genes,"Start"]
  results$End=alldist[all.genes,"End"]
  results$Strand=alldist[all.genes,"Strand"]
  results$ExonicLength=exoniclength[all.genes]
  results$NbExons=nbexons[all.genes]
  results$UnmappableLength=ovunmap[all.genes, "UnmappableLength"]

  results$HasWrongStrandJunctions=rep("No", dim(results)[1])
  results$HasWrongStrandJunctions[which(results$GeneID%in%genes.wrong.junctions)]="Yes"

  results$AllSupportedJunctions=rep("No", dim(results)[1])
  results$AllSupportedJunctions[which(results$NbExons==1)]="Yes"
  results$AllSupportedJunctions[which(results$GeneID%in%genes.all.supported)]="Yes"

  results$SomeSupportedJunctions=rep("No", dim(results)[1])
  results$SomeSupportedJunctions[which(results$NbExons==1)]="Yes"
  results$SomeSupportedJunctions[which(results$GeneID%in%genes.some.supported)]="Yes"
  
  
  results$Class=class[all.genes,"Class"]

  class.cols=c(grep("BlastX", colnames(class), value=T), grep("CSF", colnames(class), value=T))
  
  for(col in class.cols){
    results[,col]=class[all.genes, col]
  }
  
  results$Distance5PC=alldist[all.genes, "Distance5"]
  results$Distance3PC=alldist[all.genes, "Distance3"]
  results$MaxRatioSenseAntisense=maxratio.gene[results$GeneID]
  results$TotalReadCount=totreads[all.genes]
  results$ReadCountMainStrain=totreads.correctstrain[all.genes]

  results$MaxTPM=maxtpm[all.genes]

  results$ExpressionCorrelationNoOverlap=expcomp[all.genes, "ExpressionCorrelation"]
  results$ReadCountMainStrainNoOverlap=expcomp[all.genes, "ReadCountNoOverlap"]
  
  results$RegionPCSense=regions[all.genes, "RegionSense"]
  results$RegionPCAntisense=regions[all.genes, "RegionAntisense"]
  
  ## results$ProjectedPCClusters=rep("No", dim(results)[1])
  ## results$ProjectedPCClusters[which(results$GeneID%in%pc.projection.clusters$GeneID)]="Yes"

  ## results$TrinityPCClusters=rep("No", dim(results)[1])
  ## results$TrinityPCClusters[which(results$GeneID%in%pc.trinity.clusters$GeneID)]="Yes"

  results$BidirectionalPromoter1kb=rep(NA, dim(results)[1])
  results$BidirectionalPromoter1kb[which(results$GeneID%in%names(closegenes1kb))]=closegenes1kb[results$GeneID[which(results$GeneID%in%names(closegenes1kb))]]

  results$BidirectionalPromoterProteinCoding1kb=rep(NA, dim(results)[1])
  results$BidirectionalPromoterProteinCoding1kb[which(results$GeneID%in%names(closepcgenes1kb))]=closepcgenes1kb[results$GeneID[which(results$GeneID%in%names(closepcgenes1kb))]]
  
  results$BidirectionalPromoter5kb=rep(NA, dim(results)[1])
  results$BidirectionalPromoter5kb[which(results$GeneID%in%names(closegenes5kb))]=closegenes5kb[results$GeneID[which(results$GeneID%in%names(closegenes5kb))]]

  results$BidirectionalPromoterProteinCoding5kb=rep(NA, dim(results)[1])
  results$BidirectionalPromoterProteinCoding5kb[which(results$GeneID%in%names(closepcgenes5kb))]=closepcgenes5kb[results$GeneID[which(results$GeneID%in%names(closepcgenes5kb))]]

  results$CGI500bp=rep(NA, dim(results)[1])
  results$CGI500bp[which(results$GeneID%in%names(cgigenes))]=cgigenes[results$GeneID[which(results$GeneID%in%names(cgigenes))]]


  results$OverlapRNARepeats=ovrnarep[results$GeneID,"TotalFraction"]
  results$OverlapAllRepeats=ovallrep[results$GeneID,"TotalFraction"]

  results$OverlapAllRepeats[which(results$OverlapAllRepeats>1)]=1
  results$OverlapRNARepeats[which(results$OverlapRNARepeats>1)]=1

  results$OverlapRetrogenes=ovretro[results$GeneID, "TotalFraction"]

  if(length(which(results$OverlapRetrogenes>1))>0){
    stop("weird! some genes have more than 100% overlap with retrogenes")
  }

  results$miRNAPrecursor=rep("No", dim(results)[1])
  results$miRNAPrecursor[which(results$GeneID%in%mirna.precursors)]="Yes"

  results$tRNAPrecursor=rep("No", dim(results)[1])
  results$tRNAPrecursor[which(results$GeneID%in%trna.precursors)]="Yes"

  results$OverlapTSSEncodeEnhancer=oven1000.encode[all.genes]
  results$OverlapTSSVistaEnhancer=oven1000.Vista[all.genes]

  rownames(results)=results$GeneID

  for(othersp in setdiff(splist, sp)){
    ortho=get(paste("ortho.",othersp,sep=""))
    results[,paste("Ortho.",othersp,sep="")]=rep(NA, dim(results)[1])
    withortho=intersect(all.genes,  ortho[,paste("ID.",sp,sep="")])
    results[withortho, paste("Ortho.",othersp,sep="")]=ortho[withortho,paste("ID.",othersp,sep="")]
 }
  
  ## write results

  write.table(results, file=paste(pathResults, sp, "/AllInfo_FilteredTranscripts_StringTie_Ensembl", release, ".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

}

#############################################################################

