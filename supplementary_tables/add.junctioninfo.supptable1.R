######################################################################

path="LncEvoDevo/"
pathHisat=paste(path, "results/hisat/", sep="")
pathDocs=paste(path, "docs/", sep="")

######################################################################

stats.mouse=read.table(paste(pathHisat, "Mouse/stats_antisense_junctions_annotated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(stats.mouse)=paste("Mouse", stats.mouse$SampleID, sep="_")

stats.rat=read.table(paste(pathHisat, "Rat/stats_antisense_junctions_annotated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(stats.rat)=paste("Rat", stats.rat$SampleID, sep="_")

stats.chicken=read.table(paste(pathHisat, "Chicken/stats_antisense_junctions_annotated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(stats.chicken)=paste("Chicken", stats.chicken$SampleID, sep="_")

all.stats=rbind(stats.mouse, stats.rat, stats.chicken, stringsAsFactors=F)

libtype.chicken=read.table(paste(pathDocs, "LibraryType_Chicken.txt", sep=""), h=F, stringsAsFactors=F, sep="\t")

unstranded=libtype.chicken[which(libtype.chicken[,2]=="fr-unstranded"),1]

######################################################################

table=read.table("SupplementaryTable1_withoutJunctionInfo.txt", h=T, sep="\t", quote="\"", check.names=F)
table$LibraryType=rep("reverse", dim(table)[1])
table$LibraryType[which(table$Species=="Chicken" & table[,"Sample ID"]%in%unstranded)]="unstranded"

table$NbReadsCorrectSpliceJunctions=rep(0, dim(table)[1])
table$NbReadsWrongSpliceJunctions=rep(0, dim(table)[1])

table$NbReadsCorrectSpliceJunctions=all.stats[paste(table$Species, table[,"Sample ID"], sep="_"), "NbCorrect"]
table$NbReadsWrongSpliceJunctions=all.stats[paste(table$Species, table[,"Sample ID"], sep="_"), "NbWrong"]
table$NbReadsWrongSpliceJunctions[which(table$LibraryType=="unstranded")]=0

table$PropReadsWrongSpliceJunctions=100*table$NbReadsWrongSpliceJunctions/(table$NbReadsCorrectSpliceJunctions+table$NbReadsWrongSpliceJunctions)

######################################################################

colnames(table)[which(colnames(table)=="LibraryType")]="Library type"
colnames(table)[which(colnames(table)=="NbReadsCorrectSpliceJunctions")]="Nb. reads splice junctions"
colnames(table)[which(colnames(table)=="NbReadsWrongSpliceJunctions")]="Nb. spliced reads wrong strand"
colnames(table)[which(colnames(table)=="PropReadsWrongSpliceJunctions")]="Percentage spliced reads wrong strand"

write.table(table, file="SupplementaryTable1.txt", sep="\t", quote=F, row.names=F, col.names=T)

######################################################################


