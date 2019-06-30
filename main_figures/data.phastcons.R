#################################################################

source("parameters.R")
pathConservation=paste(path, "results/sequence_evolution/phastcons/Mouse/", sep="")

release=94

#################################################################

for(type in c("placental", "60vertebrates")){
  
  phastcons=list()

  exons=read.table(paste(pathConservation, type, "/PhastCons_GeneAverage_ExonBlocks_ExcludingOverlapOtherGenes_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F)
  exons$CorrectedScore=exons$Score*exons$CoveredLength/exons$UnmaskedLength
  exons[which(is.na(exons$Score)), "CorrectedScore"]=0 ## NA values means no alignment
  
  exons.score=exons$CorrectedScore
  names(exons.score)=exons$Gene

  introns=read.table(paste(pathConservation, type, "/PhastCons_IntronAverage_MaskedExons_IntronBlocks_AllGenes_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F)
  introns$CorrectedScore=introns$Score*introns$CoveredLength/introns$UnmaskedLength
  introns[which(is.na(introns$Score)), "CorrectedScore"]=0 ## NA values means no alignment
  
  introns.score=introns$CorrectedScore
  names(introns.score)=introns$Gene

  promoters=read.table(paste(pathConservation, type, "/PhastCons_Promoters_MaskedExons_400bp_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F)
  promoters$Gene=unlist(lapply(promoters$ID, function(x) unlist(strsplit(x, split="_"))[1]))
  promoters$CorrectedScore=promoters$Score*promoters$CoveredLength/promoters$UnmaskedLength
  promoters[which(is.na(promoters$Score)), "CorrectedScore"]=0
  
  promoters.score=tapply(promoters$CorrectedScore, as.factor(promoters$Gene), max, na.rm=T)
  names(promoters.score)=levels(as.factor(promoters$Gene))

  splicesites=read.table(paste(pathConservation, type, "/PhastCons_SpliceSites_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T,  stringsAsFactors=F)
  splice5.score=splicesites$SumScore5/splicesites$TotBases5
  splice3.score=splicesites$SumScore3/splicesites$TotBases3
  names(splice5.score)=splicesites$Gene
  names(splice3.score)=splicesites$Gene

  splicesites.score=(splice5.score+splice3.score)/2
  names(splicesites.score)=splicesites$Gene


  intergenic=read.table(paste(pathConservation, type, "/PhastCons_IntergenicRegions_MinDistance5kb_MinSize1kb_FilteredTranscripts_StringTie_Ensembl", release, ".txt", sep=""), h=T,  stringsAsFactors=F)
  intergenic$CorrectedScore=intergenic$Score*intergenic$CoveredLength/intergenic$UnmaskedLength
  intergenic[which(is.na(intergenic$Score)), "CorrectedScore"]=0 ## NA values means no alignment
  
  intergenic.score=intergenic$CorrectedScore
  names(intergenic.score)=intergenic$ID

  mean.intergenic=sum(intergenic$CorrectedScore*intergenic$UnmaskedLength/sum(intergenic$UnmaskedLength))
  
  phastcons[["exons"]]=exons.score
  phastcons[["introns"]]=introns.score
  phastcons[["promoters"]]=promoters.score
  phastcons[["splice5"]]=splice5.score
  phastcons[["splice3"]]=splice3.score
  phastcons[["splicesites"]]=splicesites.score
  phastcons[["intergenic"]]=intergenic.score
  phastcons[["mean.intergenic"]]=mean.intergenic

  save(phastcons, file=paste("RData/data.phastcons.Mouse.", type,".RData", sep=""))
  
}

#################################################################

