###################################################################
###################################################################

path="LncEvoDevo/"
pathProjections=paste(path, "results/exon_projections/", sep="")
pathReadCounts=paste(path, "results/expression_estimation/", sep="")
pathOrtho=paste(path, "results/ortho_genes/whole_genome_alignments/", sep="")
pathDatasets=paste(path, "supplementary_datasets/", sep="")

release=94

###################################################################
###################################################################

splist=c("Mouse", "Rat", "Chicken")

samples=list()
samples[["Mouse"]]=c("Brain_Adult1", "Brain_Adult2", "Brain_Aged1", "Brain_Aged2", "Brain_EarlyEmbryo1", "Brain_EarlyEmbryo2", "Brain_LateEmbryo1", "Brain_LateEmbryo2", "Brain_Newborn1", "Brain_Newborn2", "Kidney_Adult1", "Kidney_Adult2", "Kidney_Aged1", "Kidney_Aged2", "Kidney_Aged3", "Kidney_Aged4", "Kidney_EarlyEmbryo1", "Kidney_EarlyEmbryo2", "Kidney_LateEmbryo1", "Kidney_LateEmbryo2", "Kidney_Newborn1", "Kidney_Newborn2", "Kidney_Newborn3", "Liver_Adult1", "Liver_Adult2", "Liver_Aged1", "Liver_Aged2", "Liver_Aged3", "Liver_Aged4", "Liver_EarlyEmbryo1", "Liver_EarlyEmbryo2", "Liver_LateEmbryo1", "Liver_LateEmbryo2", "Liver_Newborn1", "Liver_Newborn2", "Liver_Newborn3", "Testis_Adult1", "Testis_Adult2", "Testis_Aged1", "Testis_Aged2", "Testis_LateEmbryo1", "Testis_LateEmbryo2", "Testis_Newborn1", "Testis_Newborn2", "Testis_Newborn3")

samples[["Rat"]]=c("Brain_Adult1", "Brain_Adult2", "Brain_Aged1", "Brain_Aged2", "Brain_EarlyEmbryo1", "Brain_EarlyEmbryo2", "Brain_LateEmbryo1", "Brain_LateEmbryo2", "Brain_Newborn1", "Brain_Newborn2", "Kidney_Adult1", "Kidney_Adult2", "Kidney_Aged1", "Kidney_Aged2", "Kidney_EarlyEmbryo1", "Kidney_EarlyEmbryo2", "Kidney_LateEmbryo1", "Kidney_LateEmbryo2", "Kidney_Newborn1", "Kidney_Newborn2", "Liver_Adult1", "Liver_Adult2", "Liver_Aged1", "Liver_Aged2", "Liver_Aged3", "Liver_Aged4", "Liver_EarlyEmbryo1", "Liver_EarlyEmbryo2", "Liver_LateEmbryo1", "Liver_LateEmbryo2", "Liver_Newborn1", "Liver_Newborn2", "Testis_Adult1", "Testis_Adult2", "Testis_Aged1", "Testis_Aged2", "Testis_LateEmbryo1", "Testis_LateEmbryo2", "Testis_Newborn1", "Testis_Newborn2")

samples[["Chicken"]]=c("Brain_EarlyEmbryo1", "Brain_EarlyEmbryo2", "Brain_LateEmbryo1", "Brain_LateEmbryo2", "Kidney_EarlyEmbryo1", "Kidney_EarlyEmbryo2", "Kidney_LateEmbryo1", "Kidney_LateEmbryo2", "Liver_EarlyEmbryo1", "Liver_EarlyEmbryo2", "Liver_LateEmbryo1", "Liver_LateEmbryo2")

###################################################################
###################################################################

for(ref in splist){
  print(ref)
  
  for(tg in setdiff(splist, ref)){
    nbsamples=length(samples[[tg]])

    print(paste(nbsamples, "samples for",tg))
    
    projstats=read.table(paste(pathProjections, "ProjectionStats_From",ref,"_To",tg,"_ExonBlocks_FilteredTranscripts_StringTie_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F)
    rownames(projstats)=projstats$GeneID

    readcounts=read.table(paste(pathReadCounts, tg, "/AllSamples_UniqueReadCounts_Projections_From",ref,"_StringTie.txt", sep=""), h=T, stringsAsFactors=F)
    rownames(readcounts)=readcounts$GeneID
    readcounts=readcounts[,-1]
    
    total.mainstrain=apply(readcounts[,samples[[tg]]], 1, sum)
    total.allsamples=apply(readcounts, 1, sum)

    names(total.mainstrain)=rownames(readcounts)
    names(total.allsamples)=rownames(readcounts)

    results=data.frame("GeneID"=projstats$GeneID, "TotalExonicLength"=projstats$TotalExonicLength, "Length.UnfilteredProjections"=projstats$UnfilteredLength, "Length.FilteredProjectionsStep1"=projstats$FilteredStep1Length,  "Length.FilteredProjectionsStep2"=projstats$FilteredStep2Length, stringsAsFactors=F)

    rownames(results)=results$GeneID

    results[,"TotalReadCount.TargetSpecies.AllSamples"]=rep(NA, dim(results)[1])
    results[rownames(readcounts), "TotalReadCount.TargetSpecies.AllSamples"]=total.allsamples
    
    results[,"TotalReadCount.TargetSpecies.MainSamples"]=rep(NA, dim(results)[1])
    results[rownames(readcounts), "TotalReadCount.TargetSpecies.MainSamples"]=total.mainstrain

    write.table(results, file=paste(pathDatasets, "SupplementaryDataset5/StatsProjections_From",ref, "_To", tg,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  }
}

###################################################################
###################################################################

alnstats.mr=read.table(paste(pathOrtho, "AlignmentStatistics_ExonBlocks_TBA_Mouse_Rat_ExonBlocks_FilteredTranscripts_StringTie_Ensembl94.txt", sep=""), h=T, stringsAsFactors=F)
rownames(alnstats.mr)=paste(alnstats.mr$ID.Mouse, alnstats.mr$ID.Rat, sep="_")

alnstats.mc=read.table(paste(pathOrtho, "AlignmentStatistics_ExonBlocks_TBA_Mouse_Chicken_ExonBlocks_FilteredTranscripts_StringTie_Ensembl94.txt", sep=""), h=T, stringsAsFactors=F)
rownames(alnstats.mc)=paste(alnstats.mc$ID.Mouse, alnstats.mc$ID.Chicken, sep="_")

alnstats.rc=read.table(paste(pathOrtho, "AlignmentStatistics_ExonBlocks_TBA_Rat_Chicken_ExonBlocks_FilteredTranscripts_StringTie_Ensembl94.txt", sep=""), h=T, stringsAsFactors=F)
rownames(alnstats.rc)=paste(alnstats.rc$ID.Rat, alnstats.rc$ID.Chicken, sep="_")


genefamilies=read.table(paste(pathDatasets, "SupplementaryDataset5/PredictedOrthoFamilies_Mouse_Rat_Chicken.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
id.mr=paste(genefamilies$ID.Mouse, genefamilies$ID.Rat, sep="_")
id.mc=paste(genefamilies$ID.Mouse, genefamilies$ID.Chicken, sep="_")
id.rc=paste(genefamilies$ID.Rat, genefamilies$ID.Chicken, sep="_")

alnstats.mr=alnstats.mr[which(rownames(alnstats.mr)%in%id.mr),]
alnstats.mc=alnstats.mc[which(rownames(alnstats.mc)%in%id.mc),]
alnstats.rc=alnstats.rc[which(rownames(alnstats.rc)%in%id.rc),]


info.mouse=read.table(paste(pathDatasets, "SupplementaryDataset1/GeneInfo_Mouse.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(info.mouse)=info.mouse$GeneID

info.rat=read.table(paste(pathDatasets, "SupplementaryDataset1/GeneInfo_Rat.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(info.rat)=info.rat$GeneID

info.chicken=read.table(paste(pathDatasets, "SupplementaryDataset1/GeneInfo_Chicken.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(info.chicken)=info.chicken$GeneID


alnstats.mr$ExonicLength.Mouse=info.mouse[alnstats.mr$ID.Mouse,"ExonicLength"]
alnstats.mr$ExonicLength.Rat=info.rat[alnstats.mr$ID.Rat,"ExonicLength"]

alnstats.mc$ExonicLength.Mouse=info.mouse[alnstats.mc$ID.Mouse,"ExonicLength"]
alnstats.mc$ExonicLength.Chicken=info.chicken[alnstats.mc$ID.Chicken,"ExonicLength"]

alnstats.rc$ExonicLength.Rat=info.rat[alnstats.rc$ID.Rat,"ExonicLength"]
alnstats.rc$ExonicLength.Chicken=info.chicken[alnstats.rc$ID.Chicken,"ExonicLength"]


write.table(alnstats.mr, file=paste(pathDatasets, "SupplementaryDataset5/AlignmentStatistics_PredictedOrthoFamilies_Mouse_Rat.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

write.table(alnstats.rc, file=paste(pathDatasets, "SupplementaryDataset5/AlignmentStatistics_PredictedOrthoFamilies_Rat_Chicken.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

write.table(alnstats.mc, file=paste(pathDatasets, "SupplementaryDataset5/AlignmentStatistics_PredictedOrthoFamilies_Mouse_Chicken.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

###################################################################
###################################################################
