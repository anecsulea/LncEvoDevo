##########################################################################

path="LncEvoDevo/"
pathLnc=paste(path, "results/lncRNA_dataset/", sep="")
pathRepeats=paste(path, "results/overlap_repeats/", sep="")

##########################################################################

sp="Mouse"
prefix="FilteredTranscripts_StringTie_Ensembl94"

##########################################################################

prom.ovrep=read.table(paste(pathRepeats, sp, "/OverlapRepeats_TEFamily_BothStrands_Promoters1kb_",prefix,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
rownames(prom.ovrep)=prom.ovrep$GeneID
prom.ovrep=prom.ovrep[,-1]
prom.ovrep$TotalFraction=apply(prom.ovrep,1,sum)
prom.ovrep$TotalFraction[which(prom.ovrep$TotalFraction>1)]=1


exons.ovrep=read.table(paste(pathRepeats, sp, "/OverlapRepeats_TEFamily_BothStrands_",prefix,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
rownames(exons.ovrep)=exons.ovrep$GeneID
exons.ovrep=exons.ovrep[,-1]
exons.ovrep$TotalFraction=apply(exons.ovrep,1,sum)
exons.ovrep$TotalFraction[which(exons.ovrep$TotalFraction>1)]=1

##########################################################################

lnc.notcons=readLines(paste(pathLnc, sp, "/SelectedLncRNAs_", prefix,"_ConservedSequence_NoExpression.txt", sep=""))
lnc.ortho=read.table(paste(pathLnc,"/PutativeOrtho_Mouse_Rat_",prefix,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
lnc.ortho=lnc.ortho$ID.Mouse

##########################################################################
