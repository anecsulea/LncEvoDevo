##############################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/",sep="")
pathTables=paste(path, "supplementary_tables/",sep="")
pathMainFigures=paste(path, "scripts/main_figures/",sep="")

##############################################################################

phast=read.table(paste(pathTables, "SupplementaryTable4.txt",sep=""), h=T, stringsAsFactors=F)
rownames(phast)=phast$GeneID

##############################################################################

expdiv=read.table(paste(pathDatasets, "SupplementaryDataset7/ExpressionDivergence_OrthoGenes_Mouse_Rat.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
expdiv$IDMouse=unlist(lapply(expdiv$ID, function(x) unlist(strsplit(x, split="_"))[1]))
expdiv$IDRat=unlist(lapply(expdiv$ID, function(x) unlist(strsplit(x, split="_"))[2]))

expdiv=expdiv[,c("IDMouse", "IDRat", "GeneType", "MeanTPM", "MaxTPM", "ExpressionDivergence")]

expdiv=expdiv[which(!is.na(expdiv$ExpressionDivergence)),]

lm1=lm(expdiv$ExpressionDivergence~log2(expdiv$MeanTPM+1))

expdiv$ResidualExpressionDivergence=lm1$residuals

##############################################################################

alnstat=read.table(paste(pathDatasets, "SupplementaryDataset5/AlignmentStatistics_PredictedOrthoFamilies_Mouse_Rat.txt", sep=""), h=T, stringsAsFactors=F)
rownames(alnstat)=alnstat$ID.Mouse

##############################################################################

allresults=expdiv
allresults$PhastConsExons=phast[allresults$IDMouse,"PhastConsExons"]
allresults$PhastConsPromoters=phast[allresults$IDMouse,"PhastConsPromoters"]
allresults$PhastConsSpliceSites=phast[allresults$IDMouse,"PhastConsSpliceSites"]

allresults$ExonicLengthMouse=alnstat[allresults$IDMouse, "ExonicLength.Mouse"]
allresults$ExonicLengthRat=alnstat[allresults$IDMouse, "ExonicLength.Rat"]
allresults$UngappedExonicLengthMouseRat=alnstat[allresults$IDMouse, "LengthUngapped"]
allresults$IdenticalExonicLengthMouseRat=alnstat[allresults$IDMouse, "LengthIdentical"]

allresults=allresults[order(allresults$GeneType),]

##############################################################################

write.table(allresults, file=paste(pathTables, "SupplementaryTable9.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

##############################################################################
