##########################################################################################

source("../main_figures/parameters.R")

##########################################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathScripts=paste(path, "scripts/expression_estimation/", sep="")

##########################################################################################

expdata=read.table(paste(pathDatasets, "SupplementaryDataset6/NormalizedTPM_OrthoGenes_Mouse_Rat.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(expdata)=expdata$ID

info=expdata[,c("ID", "GeneType")]
info$Mouse=unlist(lapply(info$ID, function(x) unlist(strsplit(x, split="_"))[1]))
info$Rat=unlist(lapply(info$ID, function(x) unlist(strsplit(x, split="_"))[2]))
  
expdata=expdata[,which(!(colnames(expdata)%in%c("ID", "GeneType")))]


##########################################################################################

avgtpm.mouse=read.table(paste(pathDatasets, "SupplementaryDataset3/Statistics_AverageTPM_Mouse.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(avgtpm.mouse)=avgtpm.mouse$GeneID
avgtpm.mouse=avgtpm.mouse[,grep("MeanTPM", colnames(avgtpm.mouse))]
colnames(avgtpm.mouse)=unlist(lapply(colnames(avgtpm.mouse), function(x) unlist(strsplit(x, split="\\."))[2]))

avgtpm.mouse=avgtpm.mouse[info$Mouse,]
rownames(avgtpm.mouse)=info$ID

##########################################################################################

avgtpm.rat=read.table(paste(pathDatasets, "SupplementaryDataset3/Statistics_AverageTPM_Rat.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(avgtpm.rat)=avgtpm.rat$GeneID
avgtpm.rat=avgtpm.rat[,grep("MeanTPM", colnames(avgtpm.rat))]
colnames(avgtpm.rat)=unlist(lapply(colnames(avgtpm.rat), function(x) unlist(strsplit(x, split="\\."))[2]))

avgtpm.rat=avgtpm.rat[info$Rat,]
rownames(avgtpm.rat)=info$ID

avgtpm.rat=avgtpm.rat[,colnames(avgtpm.mouse)]

##########################################################################################

minexp=1
minsamples=2

all.genes=info$ID

##########################################################################################

results=list()

samples=colnames(avgtpm.mouse)

allexpdata=expdata
expdata.mouse=avgtpm.mouse
expdata.rat=avgtpm.rat

nbsamplesexp.mouse=apply(allexpdata[,grep("Mouse", colnames(allexpdata))], 1, function(x) length(which(x>=minexp)))
nbsamplesexp.rat=apply(allexpdata[,grep("Rat", colnames(allexpdata))], 1, function(x) length(which(x>=minexp)))
selected.genes=rownames(allexpdata)[which(nbsamplesexp.mouse>=minsamples & nbsamplesexp.rat>=minsamples)]

## take only selected genes

expdata.mouse=expdata.mouse[selected.genes,]
expdata.rat=expdata.rat[selected.genes,]

## divide by sum expression

relexpdata.mouse=expdata.mouse/apply(expdata.mouse,1,sum)
relexpdata.rat=expdata.rat/apply(expdata.rat,1,sum)

## compute distance

distance=unlist(lapply(1:dim(expdata.mouse)[1], function(x) {y=as.numeric(relexpdata.mouse[x,]); z=as.numeric(relexpdata.rat[x,]); return(sqrt(sum((y-z)^2)))}))
names(distance)=selected.genes

contributions=lapply(samples, function(x) ((relexpdata.mouse[,x]-relexpdata.rat[,x])^2)/(distance^2))
names(contributions)=samples


meanexp=apply(cbind(expdata.mouse, expdata.rat),1, mean)
maxexp=apply(cbind(expdata.mouse, expdata.rat),1, max)

results[["MeanTPM"]]=meanexp
results[["MaxTPM"]]=maxexp
results[["ExpressionDivergence"]]=distance

for(sample in samples){
  results[[paste("PCExpDiv.",sample,sep="")]]=100*contributions[[sample]]
}


results=as.data.frame(results)
results$ID=selected.genes
results$GeneType=info[selected.genes,"GeneType"]

results=results[,c("ID", "GeneType", setdiff(colnames(results), c("ID", "GeneType")))]


##########################################################################################

write.table(results, file=paste(pathDatasets, "SupplementaryDataset8/ExpressionDivergence_OrthoGenes_Mouse_Rat.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

##########################################################################################
