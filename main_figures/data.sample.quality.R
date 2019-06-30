###############################################################################

source("parameters.R")

pathHisat=paste(path, "results/hisat/", sep="")
pathDegradation=paste(path, "results/RNA_degradation/", sep="")

###############################################################################

deg=read.table(paste(pathDegradation, "Stats_RIN_3bias.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

deg$Sample.ID=paste(deg$Species, deg$Sample.ID, sep="_")

rownames(deg)=deg$Sample.ID

colnames(deg)=c("Species", "SampleID", "NbReads", "RIN", "CoverageBias")

###############################################################################

save(deg, file="RData/data.sample.quality.RData")

###############################################################################
