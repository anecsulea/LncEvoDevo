##################################################################################

source("parameters.R")

pathEnsembl=paste(path, "data/ensembl_annotations/", sep="")

release=94

##################################################################################

for(sp in c("Mouse", "Rat", "Chicken")){
  print(sp)
  
  go=read.table(paste(pathEnsembl, sp, "/GO_Ensembl94.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

  print(paste(dim(go)[1], "lines in the file"))

  go=go[which(go$GO.term.accession!=""),]

  go=go[which(!(go$GO.term.name%in%c("biological_process", "cellular_component", "molecular_function"))),] ## remove the generic processes

  GOdata=list()

  for(domain in c("biological_process", "cellular_component", "molecular_function")){
    GOdata[[domain]]=list()
    
    this.domain=go[which(go$GO.domain==domain),]
    categories=this.domain[,c(2,3)]

    dupli=which(duplicated(categories[,1]))
    if(length(dupli)>0){
      categories=categories[-dupli,]
    }

    rownames(categories)=categories[,1]

    genecat=this.domain[,c(1,2)]

    genelist=lapply(unique(genecat$GO.term.accession), function(x) genecat[which(genecat[,2]==x),1])
    names(genelist)=unique(genecat$GO.term.accession)

    GOdata[[domain]][["categories"]]=categories
    GOdata[[domain]][["genelist"]]=genelist
  }
  
  save(GOdata, file=paste("RData/data.gene.ontology.",sp,".RData", sep=""))
}

##################################################################################


