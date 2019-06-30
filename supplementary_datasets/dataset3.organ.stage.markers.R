##########################################################################

path="LncEvoDevo/"
pathExpression=paste(path, "results/expression_estimation/", sep="")
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathSupplementaryTables=paste(path, "supplementary_tables/", sep="")
pathMainFigures=paste(path, "scripts/main_figures/", sep="")

options(scipen=999) ## remove scientific notation ## options(scipen=0) to get it back

release=94

minTPM=2
mintau=0.85
maxFDRGO=0.1

###############################################################################

source(paste(pathMainFigures, "compute.go.enrichment.R",sep=""))

load(paste(pathMainFigures, "RData/data.gene.ontology.Mouse.RData",sep=""))

###############################################################################

info.mouse=read.table(paste(pathDatasets, "SupplementaryDataset1/GeneInfo_Mouse.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(info.mouse)=info.mouse$GeneID

###############################################################################

stats.mouse=read.table(paste(pathDatasets, "SupplementaryDataset3/Statistics_AverageTPM_Mouse.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(stats.mouse)=stats.mouse$GeneID

stats.rat=read.table(paste(pathDatasets, "SupplementaryDataset3/Statistics_AverageTPM_Rat.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(stats.rat)=stats.rat$GeneID

###############################################################################

ortho.mr=read.table(paste(pathDatasets, "SupplementaryDataset5/EnsemblOrtho_ProteinCodingGenes_1to1_Mouse_Rat.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
ortho.mr=ortho.mr[which(ortho.mr$Mouse%in%stats.mouse$GeneID & ortho.mr$Rat%in%stats.rat$GeneID),]

###############################################################################

common.expdata=ortho.mr[,c("Mouse", "Rat")]
colnames(common.expdata)=c("ID.Mouse", "ID.Rat")

common.expdata$MaxSample.Mouse=stats.mouse[common.expdata$ID.Mouse, "MaxSample4Stages"]
common.expdata$MaxExpression.Mouse=stats.mouse[common.expdata$ID.Mouse, "MaxExpression4Stages"]
common.expdata$ExpressionSpecificity.Mouse=stats.mouse[common.expdata$ID.Mouse, "ExpressionSpecificity4Stages"]

common.expdata$MaxSample.Rat=stats.rat[common.expdata$ID.Rat, "MaxSample4Stages"]
common.expdata$MaxExpression.Rat=stats.rat[common.expdata$ID.Rat, "MaxExpression4Stages"]
common.expdata$ExpressionSpecificity.Rat=stats.rat[common.expdata$ID.Rat, "ExpressionSpecificity4Stages"]

###############################################################################

filtered.expdata=common.expdata[which(common.expdata$MaxExpression.Mouse>=minTPM & common.expdata$MaxExpression.Rat>=minTPM & common.expdata$ExpressionSpecificity.Mouse>=mintau & common.expdata$ExpressionSpecificity.Rat>=mintau),]
filtered.expdata=filtered.expdata[order(filtered.expdata$ExpressionSpecificity.Rat, decreasing=T),]
filtered.expdata=filtered.expdata[order(filtered.expdata$ExpressionSpecificity.Mouse, decreasing=T),]

nbshared=length(which(filtered.expdata$MaxSample.Mouse==filtered.expdata$MaxSample.Rat))

nbtot=dim(filtered.expdata)[1]

print(paste(nbshared,"shared specificity out of", nbtot,  round(100*nbshared/nbtot),"%"))

shared=filtered.expdata[which(filtered.expdata$MaxSample.Mouse==filtered.expdata$MaxSample.Rat),]
shared$GeneName=info.mouse[shared$ID.Mouse, "GeneName"]
shared=shared[,c(c("ID.Mouse", "ID.Rat", "GeneName"), setdiff(colnames(shared), c("ID.Mouse", "ID.Rat", "GeneName")))]
## this is also supplementary table 4!

write.table(shared, paste(pathSupplementaryTables, "SupplementaryTable4.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
            
###############################################################################

backgroundlist=ortho.mr$Mouse

for(sample in unique(shared$MaxSample.Mouse)){
  print(sample)
  
  this.expdata=shared[which(shared$MaxSample.Mouse==sample),]
  write.table(this.expdata, file=paste(pathDatasets, "SupplementaryDataset3/GeneList_OrganStageMarkers_MouseRat_",sample,".txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

  this.go=compute.go.enrichment(this.expdata$ID.Mouse, backgroundlist, GOdata[["biological_process"]][["categories"]], GOdata[["biological_process"]][["genelist"]])
  this.go=this.go[order(this.go$FDR),]
  this.go=this.go[which(this.go$FDR<maxFDRGO),]

  if(dim(this.go)[1]>0){
    print(head(this.go[,c("GOName", "Enrichment","FDR")]))
    write.table(this.go, file=paste(pathDatasets, "SupplementaryDataset3/GOEnrichment_BiologicalProcess_OrganStageMarkers_MouseRat_",sample,".txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  }


}

###############################################################################
