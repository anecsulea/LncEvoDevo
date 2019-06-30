##############################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/",sep="")
pathTables=paste(path, "supplementary_tables/",sep="")
pathMainFigures=paste(path, "scripts/main_figures/",sep="")

##############################################################################

ortho=read.table(paste(pathDatasets,"SupplementaryDataset5/PredictedOrthoFamilies_Mouse_Rat_Chicken.txt", sep=""), h=T,stringsAsFactors=F, sep="\t", quote="\"")


mrc=ortho[which(ortho$GeneType.Mouse=="lncRNA" & ortho$GeneType.Rat=="lncRNA" & ortho$GeneType.Chicken=="lncRNA"),]

mrc=mrc[,c("ID.Mouse", "ID.Rat", "ID.Chicken")]

write.table(mrc, file=paste(pathTables, "SupplementaryTable7.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

##############################################################################
