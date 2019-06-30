##########################################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathScripts=paste(path, "scripts/expression_estimation/", sep="")

##########################################################################################

source(paste(pathScripts, "normalization.R", sep=""))

##########################################################################################

tpm.mouse=read.table(paste(pathDatasets, "SupplementaryDataset2/KallistoRawPM_Mouse.txt", sep=""), h=T, sep="\t", stringsAsFactors=F)
rownames(tpm.mouse)=tpm.mouse$GeneID
tpm.mouse=tpm.mouse[,-1]
colnames(tpm.mouse)=paste("Mouse", colnames(tpm.mouse), sep="_")

tpm.rat=read.table(paste(pathDatasets, "SupplementaryDataset2/KallistoRawPM_Rat.txt", sep=""), h=T, sep="\t", stringsAsFactors=F)
rownames(tpm.rat)=tpm.rat$GeneID
tpm.rat=tpm.rat[,-1]
colnames(tpm.rat)=paste("Rat", colnames(tpm.rat), sep="_")


tpm.chicken=read.table(paste(pathDatasets, "SupplementaryDataset2/KallistoRawPM_Chicken.txt", sep=""), h=T, sep="\t", stringsAsFactors=F)
rownames(tpm.chicken)=tpm.chicken$GeneID
tpm.chicken=tpm.chicken[,-1]
colnames(tpm.chicken)=paste("Chicken", colnames(tpm.chicken), sep="_")

##########################################################################################

## predicted ortho with our method

predicted.ortho=read.table(paste(pathDatasets, "SupplementaryDataset5/PredictedOrthoFamilies_Mouse_Rat_Chicken.txt", sep=""), h=T, sep="\t", stringsAsFactors=F)

########################################################################################

lncortho.mr=predicted.ortho[which(predicted.ortho$GeneType.Mouse=="lncRNA" & predicted.ortho$GeneType.Rat=="lncRNA"),c("ID.Mouse", "ID.Rat")]
colnames(lncortho.mr)=c("Mouse", "Rat")

##########################################################################################

lncortho.mrc=predicted.ortho[which(predicted.ortho$GeneType.Mouse=="lncRNA" & predicted.ortho$GeneType.Rat=="lncRNA"  & predicted.ortho$GeneType.Chicken=="lncRNA"),c("ID.Mouse", "ID.Rat", "ID.Chicken")]
colnames(lncortho.mrc)=c("Mouse", "Rat", "Chicken")

##########################################################################################

predicted.mr=apply(predicted.ortho[,c("ID.Mouse", "ID.Rat")],1, function(x) paste(x, collapse="_"))
predicted.mrc=apply(predicted.ortho[,c("ID.Mouse", "ID.Rat", "ID.Chicken")],1, function(x) paste(x, collapse="_"))

##########################################################################################

pcortho.mr=read.table(paste(pathDatasets, "SupplementaryDataset5/EnsemblOrtho_ProteinCodingGenes_1to1_Mouse_Rat.txt", sep=""), h=T, sep="\t", stringsAsFactors=F)
pcortho.mrc=read.table(paste(pathDatasets, "SupplementaryDataset5/EnsemblOrtho_ProteinCodingGenes_1to1_Mouse_Rat_Chicken.txt", sep=""), h=T, sep="\t", stringsAsFactors=F)

idpc.mr=apply(pcortho.mr[,c("Mouse", "Rat")],1, function(x) paste(x, collapse="_"))
idpc.mrc=apply(pcortho.mrc[,c("Mouse", "Rat", "Chicken")],1, function(x) paste(x, collapse="_"))

print(paste(length(idpc.mr),"mouse-rat ortho in Ensembl",length(intersect(predicted.mr, idpc.mr)), "confirmed"))
print(paste(length(idpc.mrc),"mouse-rat-chicken ortho in Ensembl",length(intersect(predicted.mrc, idpc.mrc)), "confirmed"))

## we take only confirmed ortho, both Ensembl and our predictions

pcortho.mr=pcortho.mr[which(idpc.mr%in%predicted.mr),]
pcortho.mrc=pcortho.mrc[which(idpc.mrc%in%predicted.mrc),]

##########################################################################################

allortho.mr=rbind(pcortho.mr[,c("Mouse", "Rat")], lncortho.mr[,c("Mouse", "Rat")], stringsAsFactors=F)
genetype.mr=c(rep("protein_coding", dim(pcortho.mr)[1]), rep("lncRNA", dim(lncortho.mr)[1]))

allortho.mrc=rbind(pcortho.mrc[,c("Mouse", "Rat", "Chicken")], lncortho.mrc[,c("Mouse", "Rat", "Chicken")], stringsAsFactors=F)
genetype.mrc=c(rep("protein_coding", dim(pcortho.mrc)[1]), rep("lncRNA", dim(lncortho.mrc)[1]))

##########################################################################################

## combine tpm values

tpm.mr=data.frame(tpm.mouse[allortho.mr$Mouse,], tpm.rat[allortho.mr$Rat,], stringsAsFactors=F)
rownames(tpm.mr)=paste(allortho.mr$Mouse, allortho.mr$Rat, sep="_")

tpm.mrc=data.frame(tpm.mouse[allortho.mrc$Mouse,], tpm.rat[allortho.mrc$Rat,],tpm.chicken[allortho.mrc$Chicken,], stringsAsFactors=F)
rownames(tpm.mrc)=paste(allortho.mrc$Mouse, allortho.mrc$Rat, allortho.mrc$Chicken, sep="_")

##########################################################################################

## normalize by housekeeping genes

normdata.mr=normalization(tpm.mr)
normdata.mrc=normalization(tpm.mrc)

normtpm.mr=normdata.mr[["expdata.norm"]]
normtpm.mrc=normdata.mrc[["expdata.norm"]]

##########################################################################################

tpm.mr$ID=rownames(tpm.mr)
tpm.mr$GeneType=genetype.mr
tpm.mr=tpm.mr[,c("ID", "GeneType", setdiff(colnames(tpm.mr), c("ID", "GeneType")))]

tpm.mrc$ID=rownames(tpm.mrc)
tpm.mrc$GeneType=genetype.mrc
tpm.mrc=tpm.mrc[,c("ID", "GeneType", setdiff(colnames(tpm.mrc), c("ID", "GeneType")))]

##########################################################################################

normtpm.mr$ID=rownames(tpm.mr)
normtpm.mr$GeneType=genetype.mr
normtpm.mr=normtpm.mr[,c("ID", "GeneType", setdiff(colnames(normtpm.mr), c("ID", "GeneType")))]

normtpm.mrc$ID=rownames(tpm.mrc)
normtpm.mrc$GeneType=genetype.mrc
normtpm.mrc=normtpm.mrc[,c("ID", "GeneType", setdiff(colnames(normtpm.mrc), c("ID", "GeneType")))]

##########################################################################################

write.table(tpm.mr, file=paste(pathDatasets, "SupplementaryDataset6/RawTPM_OrthoGenes_Mouse_Rat.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
write.table(tpm.mrc, file=paste(pathDatasets, "SupplementaryDataset6/RawTPM_OrthoGenes_Mouse_Rat_Chicken.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

write.table(normtpm.mr, file=paste(pathDatasets, "SupplementaryDataset6/NormalizedTPM_OrthoGenes_Mouse_Rat.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
write.table(normtpm.mrc, file=paste(pathDatasets, "SupplementaryDataset6/NormalizedTPM_OrthoGenes_Mouse_Rat_Chicken.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

##########################################################################################
