################################################################################

path="LncEvoDevo/"
pathSingleCell=paste(path, "data/single_cell_expression/", sep="")
pathEnsemblOrtho=paste(path, "data/ensembl_ortho/",sep="")
pathDatasets=paste(path, "supplementary_datasets/",sep="")
pathTables=paste(path, "supplementary_tables/",sep="")

################################################################################

geneinfo=read.table(paste(pathDatasets, "SupplementaryDataset1/GeneInfo_Mouse.txt", sep=""),  h=T, stringsAsFactors=F, sep="\t", quote="\"")

################################################################################

markers=read.table(paste(pathSingleCell, "TabulaMuris_CellMarkers.txt", sep=""), h=F, stringsAsFactors=F, sep="\t")

################################################################################

all.tissues=c()
all.cells=c()
all.ids=c()
all.names=c()

for(tissue in unique(markers$V1)){
  cells=unique(markers$V2[which(markers$V1==tissue)])

  for(cell in cells){

    names=unlist(strsplit(markers$V3[which(markers$V1==tissue & markers$V2==cell)], split=","))
    names=intersect(names, geneinfo$GeneName)
    
    ids=unlist(lapply(names, function(x) {y=which(geneinfo$GeneName==x); if(length(y)==1){return(geneinfo$GeneID[y])} else{return(NA)}}))
    
    all.tissues=c(all.tissues, rep(tissue, length(names)))
    all.cells=c(all.cells, rep(cell, length(names)))
    all.names=c(all.names, names)
    all.ids=c(all.ids, ids)
    
  }
}

markers=data.frame("Organ"=all.tissues, "CellType"=all.cells, "GeneName"=all.names, "MouseGeneID"=all.ids, stringsAsFactors=F)
markers=markers[which(markers$CellType!="unknown"),]
markers=markers[which(markers$CellType!="kidney_cell"),] ## we don't know what that is

dupligenes=markers$MouseGeneID[which(duplicated(markers$MouseGeneID))]

if(length(dupligenes)>0){ ## we remove genes that are said to be markers for several different cell types
  markers=markers[-which(markers$MouseGeneID%in%dupligenes),]
}

markers$CellType[grep("^kidney", markers$CellType)]=unlist(lapply(markers$CellType[grep("^kidney", markers$CellType)], function(x) substr(x, 8, nchar(x))))

markers$Reference=rep("TabulaMuris, Nature, 2018",dim(markers)[1])
markers$Reference[which(markers$Organ=="Testis")]="Green et al., Dev Cell, 2018"

markers$ImmuneCells=rep("no", dim(markers)[1])
markers$ImmuneCells[which(markers$CellType%in%c("macrophage", "leukocyte", "B_cell","Kupffer_cell", "natural_killer_cell", "innate_lymphoid"))]="yes"

################################################################################

## now add ortho genes for mouse-rat and mouse-chicken

ortho.mr=read.table(paste(pathEnsemblOrtho, "GeneFamilies_1to1_MouseRat_ProteinCoding_Ensembl94.txt", sep=""), h=T, stringsAsFactors=F)
rownames(ortho.mr)=ortho.mr$Mouse

ortho.mc=read.table(paste(pathEnsemblOrtho, "GeneFamilies_1to1_MouseChicken_ProteinCoding_Ensembl94.txt", sep=""), h=T, stringsAsFactors=F)
rownames(ortho.mc)=ortho.mc$Mouse

################################################################################

markers$RatGeneID=rep(NA, dim(markers)[1])
markers[which(markers$MouseGeneID%in%ortho.mr$Mouse),"RatGeneID"]=ortho.mr[markers$MouseGeneID[which(markers$MouseGeneID%in%ortho.mr$Mouse)], "Rat"]

markers$ChickenGeneID=rep(NA, dim(markers)[1])
markers[which(markers$MouseGeneID%in%ortho.mc$Mouse),"ChickenGeneID"]=ortho.mc[markers$MouseGeneID[which(markers$MouseGeneID%in%ortho.mc$Mouse)], "Chicken"]

################################################################################

write.table(markers, file=paste(pathTables, "SupplementaryTable3.txt", sep=""), row.names=F,col.names=T, sep="\t", quote=F)

################################################################################
