##############################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/",sep="")
pathTables=paste(path, "supplementary_tables/",sep="")
pathMainFigures=paste(path, "scripts/main_figures/",sep="")

##############################################################################

expdiv=read.table(paste(pathDatasets, "SupplementaryDataset8/ExpressionDivergence_OrthoGenes_Mouse_Rat.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
expdiv$IDMouse=unlist(lapply(expdiv$ID, function(x) unlist(strsplit(x, split="_"))[1]))
expdiv$IDRat=unlist(lapply(expdiv$ID, function(x) unlist(strsplit(x, split="_"))[2]))

expdiv=expdiv[,c("IDMouse", "IDRat", "GeneType", "MeanTPM", "MaxTPM", "ExpressionDivergence")]

expdiv=expdiv[which(!is.na(expdiv$ExpressionDivergence)),]

lm1=lm(expdiv$ExpressionDivergence~log2(expdiv$MeanTPM+1))

expdiv$ResidualExpressionDivergence=lm1$residuals

##############################################################################

expdiv=expdiv[order(expdiv$ResidualExpressionDivergence, decreasing=T),]
expdiv.pc=expdiv[which(expdiv$GeneType=="protein_coding"),]

nbtot=dim(expdiv.pc)[1]
top10pc=expdiv.pc$IDMouse[1:floor(nbtot/10)]
allpc=expdiv.pc$IDMouse

##############################################################################

source(paste(pathMainFigures, "compute.go.enrichment.R",sep=""))

load(paste(pathMainFigures,"RData/data.gene.ontology.Mouse.RData", sep=""))

this.go=compute.go.enrichment(top10pc, allpc, GOdata[["biological_process"]][["categories"]], GOdata[["biological_process"]][["genelist"]])
this.go=this.go[order(this.go$FDR),]

this.go=this.go[which(this.go$FDR<0.1),]

write.table(this.go, file=paste(pathDatasets, "SupplementaryDataset8/GOEnrichment_BiologicalProcess_Top10Percent_ProteinCodingGenes.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

##############################################################################
