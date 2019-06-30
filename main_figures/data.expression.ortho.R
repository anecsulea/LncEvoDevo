########################################################################

source("parameters.R")

pathDatasets=paste(path, "supplementary_datasets/", sep="")

########################################################################

exportho.mr=read.table(paste(pathDatasets, "SupplementaryDataset6/NormalizedTPM_OrthoGenes_Mouse_Rat.txt", sep=""), h=T, stringsAsFactors=F)
exportho.mrc=read.table(paste(pathDatasets, "SupplementaryDataset6/NormalizedTPM_OrthoGenes_Mouse_Rat_Chicken.txt", sep=""), h=T, stringsAsFactors=F)

rownames(exportho.mr)=exportho.mr$ID
rownames(exportho.mrc)=exportho.mrc$ID

########################################################################

info.mr=exportho.mr[,c("ID", "GeneType")]
info.mrc=exportho.mrc[,c("ID", "GeneType")]

tpm.mr=exportho.mr[,-c(1,2)]
tpm.mrc=exportho.mrc[,-c(1,2)]

sptissage.mr=as.factor(unlist(lapply(colnames(tpm.mr), function(x) substr(x,1, nchar(x)-1))))
sptissage.mrc=as.factor(unlist(lapply(colnames(tpm.mrc), function(x) substr(x,1, nchar(x)-1))))

avgexp.mr=t(apply(tpm.mr, 1, function(x) tapply(as.numeric(x), sptissage.mr, mean)))
avgexp.mrc=t(apply(tpm.mrc, 1, function(x) tapply(as.numeric(x), sptissage.mrc, mean)))

avgexp.mr=data.frame(info.mr, avgexp.mr)
avgexp.mrc=data.frame(info.mrc, avgexp.mrc)

########################################################################

save(list=c("exportho.mr", "exportho.mrc", "avgexp.mr", "avgexp.mrc"), file="RData/data.expression.ortho.RData")

########################################################################


