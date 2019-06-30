################################################################################

source("parameters.R")

pathDatasets=paste(path, "supplementary_datasets/", sep="")

################################################################################ 

minreads=100

mouse.specific=read.table(paste(pathDatasets, "SupplementaryDataset8/Candidate_MouseSpecific_LncRNAs_Min",minreads,"Reads.txt", sep=""), h=T, stringsAsFactors=F,sep="\t")$GeneID
rat.specific=read.table(paste(pathDatasets, "SupplementaryDataset8/Candidate_RatSpecific_LncRNAs_Min",minreads,"Reads.txt", sep=""), h=T, stringsAsFactors=F,sep="\t")$GeneID

################################################################################

predicted.ortho=read.table(paste(pathDatasets, "SupplementaryDataset5/PredictedOrthoFamilies_Mouse_Rat_Chicken.txt", sep=""), h=T, sep="\t", stringsAsFactors=F)


mr.ortho=predicted.ortho[which(predicted.ortho$GeneType.Mouse=="lncRNA" & predicted.ortho$GeneType.Rat=="lncRNA"),c("ID.Mouse", "ID.Rat")]
colnames(mr.ortho)=c("Mouse", "Rat")

mouse.ortho=mr.ortho$Mouse
rat.ortho=mr.ortho$Rat

################################################################################

load("RData/data.annotations.Mouse.RData")
allinfo.mouse=allinfo
pc.mouse=pc
lnc.mouse=lnc


load("RData/data.annotations.Rat.RData")
allinfo.rat=allinfo
pc.rat=pc
lnc.rat=lnc

mouse.new=allinfo.mouse$GeneID[which(allinfo.mouse$AnnotationSource=="DeNovo")]
mouse.known=allinfo.mouse$GeneID[which(allinfo.mouse$AnnotationSource=="Ensembl")]

rat.new=allinfo.rat$GeneID[which(allinfo.rat$AnnotationSource=="DeNovo")]
rat.known=allinfo.rat$GeneID[which(allinfo.rat$AnnotationSource=="Ensembl")]

################################################################################

save(list=c("mouse.new", "mouse.known", "rat.new", "rat.known","mouse.specific", "rat.specific", "mouse.ortho", "rat.ortho"), file="RData/data.comparison.specific.ortho.RData")

################################################################################

