###################################################################################

source("parameters.R")
pathDatasets=paste(path, "supplementary_datasets/", sep="")

###################################################################################

source("parameters.R")
source("compute.go.enrichment.R")

###################################################################################

maxFDR=0.01
maxFDRGO=0.1

signif=list()
GO=list()

for(sp in c("Mouse", "Rat")){
  signif[[sp]]=list()

  load(paste("RData/data.diffexp.", sp, ".RData", sep=""))
  de.global=de.global.allreads

  for(tiss in tissue.order){
    signif[[sp]][[tiss]]=de.global$GeneID[which(de.global[,paste("FDR", tiss, sep=".")]<maxFDR)]
  }
  
  load(paste("RData/data.gene.ontology.",sp,".RData", sep=""))

  GO[[sp]]=GOdata
}

###################################################################################

load("RData/data.expression.ortho.RData")
avgexp.mr=avgexp.mr[which(avgexp.mr$GeneType=="protein_coding"),]
avgexp.mr$IDMouse=unlist(lapply(avgexp.mr$ID, function(x) unlist(strsplit(x, split="_"))[1]))
avgexp.mr$IDRat=unlist(lapply(avgexp.mr$ID, function(x) unlist(strsplit(x, split="_"))[2]))

###################################################################################

kmeans.results=list()

for(tiss in tissue.order){
  nbclust=5

  if(tiss=="Testis"){
    nbclust=4
  }
  
  selected.genes=avgexp.mr$ID[which(avgexp.mr$IDMouse%in%signif[["Mouse"]][[tiss]] & avgexp.mr$IDRat%in%signif[["Rat"]][[tiss]])]
  
  exp.mouse=avgexp.mr[selected.genes,grep(paste("Mouse",tiss, sep="_"), colnames(avgexp.mr))]
  exp.rat=avgexp.mr[selected.genes,grep(paste("Rat", tiss, sep="_"), colnames(avgexp.mr))]
    
  exp.mouse=exp.mouse/apply(exp.mouse,1, max)
  exp.rat=exp.rat/apply(exp.rat,1, max)
  
  exp.ortho=cbind(exp.mouse, exp.rat, stringsAsFactors=F)

  has.na=apply(exp.ortho,1, function(x) any(is.na(x)))
  exp.ortho=exp.ortho[which(!has.na),]

  this.age.order=age.order
  if(tiss=="Testis"){
    this.age.order=setdiff(age.order, "EarlyEmbryo")
  }
  
  sample.order=kronecker(paste(tiss, this.age.order, sep="_"),c("Mouse", "Rat"),  function(x,y) paste(y,x, sep="_"))
  
  exp.ortho=exp.ortho[,sample.order]
  
  this.kmeans=kmeans(exp.ortho, centers=nbclust, iter.max=50)
  kmeans.results[[tiss]]=this.kmeans

  ## write data for supplementary dataset
  
  all.results=exp.ortho
  all.results$IDMouse=unlist(lapply(rownames(exp.ortho), function(x) unlist(strsplit(x, split="_"))[1]))
  all.results$IDRat=unlist(lapply(rownames(exp.ortho), function(x) unlist(strsplit(x, split="_"))[2]))
  all.results$ClusterIndex=this.kmeans$cluster
  all.results=all.results[,c(c("ClusterIndex", "IDMouse", "IDRat"), setdiff(colnames(all.results), c("ClusterIndex", "IDMouse", "IDRat")))]  

  all.genes=all.results$IDMouse
  
  for(i in unique(all.results$ClusterIndex)){
    targetlist=all.results$IDMouse[which(all.results$ClusterIndex==i)]
    this.go=compute.go.enrichment(targetlist, all.genes, GO[["Mouse"]][["biological_process"]][["categories"]], GO[["Mouse"]][["biological_process"]][["genelist"]])
    this.go=this.go[order(this.go$FDR),]


    if(length(which(this.go$FDR<maxFDRGO))>0){
      this.go=this.go[which(this.go$FDR<maxFDRGO),]

      write.table(this.go, file=paste(pathDatasets, "SupplementaryDataset4/GOEnrichment_BiologicalProcess_KmeansClusters_MouseRatOrtho_DiffExp_",tiss,"_Cluster", i, "_ProteinCoding.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
    }
  }

  write.table(all.results, file=paste(pathDatasets, "SupplementaryDataset4/KmeansClusters_MouseRatOrtho_DiffExp_",tiss, "_ProteinCoding.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

}

save(kmeans.results, file="RData/data.kmeans.diffexp.orthopc.RData")

###################################################################################
