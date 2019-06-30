########################################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/", sep="")

########################################################################################

info.mouse=read.table(paste(pathDatasets, "SupplementaryDataset1/GeneInfo_Mouse.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(info.mouse)=info.mouse$GeneID

info.rat=read.table(paste(pathDatasets, "SupplementaryDataset1/GeneInfo_Rat.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(info.rat)=info.rat$GeneID

########################################################################################

proj.mouse.rat=read.table(paste(pathDatasets, "SupplementaryDataset5/StatsProjections_FromMouse_ToRat.txt", sep=""), h=T, stringsAsFactors=F)
rownames(proj.mouse.rat)=proj.mouse.rat$GeneID
proj.mouse.rat=proj.mouse.rat[info.mouse$GeneID,]

proj.rat.mouse=read.table(paste(pathDatasets, "SupplementaryDataset5/StatsProjections_FromRat_ToMouse.txt", sep=""), h=T, stringsAsFactors=F)
rownames(proj.rat.mouse)=proj.rat.mouse$GeneID
proj.rat.mouse=proj.rat.mouse[info.rat$GeneID,]

########################################################################################

families=read.table(paste(pathDatasets, "SupplementaryDataset5/PredictedOrthoFamilies_Mouse_Rat_Chicken.txt", sep=""), h=T, stringsAsFactors=F)

########################################################################################

for(minreads in c(100, 250)){
  mouse.specific=info.mouse$GeneID[which(info.mouse$SelectedLncRNA=="Yes" &  (!info.mouse$GeneID%in%families$ID.Mouse) & info.mouse$ReadCount>=minreads & proj.mouse.rat$Length.FilteredProjectionsStep2==proj.mouse.rat$TotalExonicLength & proj.mouse.rat$TotalReadCount.TargetSpecies.AllSamples==0)]
  
  rat.specific=info.rat$GeneID[which(info.rat$SelectedLncRNA=="Yes" &  (!info.rat$GeneID%in%families$ID.Rat) & info.rat$ReadCount>=minreads & proj.rat.mouse$Length.FilteredProjectionsStep2==proj.rat.mouse$TotalExonicLength & proj.rat.mouse$TotalReadCount.TargetSpecies.AllSamples==0)]
  
  write.table(info.mouse[mouse.specific, ], file=paste(pathDatasets, "SupplementaryDataset9/Candidate_MouseSpecific_LncRNAs_Min",minreads,"Reads.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  
  write.table(info.rat[rat.specific, ], file=paste(pathDatasets, "SupplementaryDataset9/Candidate_RatSpecific_LncRNAs_Min",minreads,"Reads.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
}

########################################################################################
