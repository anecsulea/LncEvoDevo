##############################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/",sep="")
pathTables=paste(path, "supplementary_tables/",sep="")
pathMainFigures=paste(path, "scripts/main_figures/",sep="")

##############################################################################

info=read.table(paste(pathDatasets, "SupplementaryDataset1/GeneInfo_Mouse.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
pc=info$GeneID[which(info$SelectedProteinCoding=="Yes")]
lnc=info$GeneID[which(info$CandidateLncRNA=="Yes")]

load(paste(pathMainFigures,"RData/data.phastcons.Mouse.placental.RData",sep=""))

##############################################################################

genes=c(pc, lnc)
genetype=rep(c("protein_coding", "lncRNA"), c(length(pc), length(lnc)))
prom=phastcons[["promoters"]][genes]
exons=phastcons[["exons"]][genes]
introns=phastcons[["introns"]][genes]
splicesites=phastcons[["splicesites"]][genes]


results=data.frame("GeneID"=genes, "GeneType"=genetype,"PhastConsExons"=exons, "PhastConsPromoters"=prom, "PhastConsSpliceSites"=splicesites, "PhastConsIntrons"=introns, stringsAsFactors=F)

write.table(results, file=paste(pathTables, "SupplementaryTable6.txt", sep=""),row.names=F, col.names=T, sep="\t", quote=F)

##############################################################################
