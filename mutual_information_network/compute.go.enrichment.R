################################################################################

path="/sps/biometr/necsulea/LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathResults=paste(path, "results/mutual_information_network/", sep="")
pathFigures=paste(path, "scripts/main_figures/",sep="")

################################################################################

options(scipen=999)

maxFDR=0.001 ## to define mutual information network

type="lncRNAs_only"

minConnections=50 ## minimum number of connected protein-coding genes to compute GO enrichment

maxGOFDR=0.1

############################################################################

source(paste(pathFigures, "compute.go.enrichment.R", sep=""))

load(paste(pathFigures, "RData/data.gene.ontology.Mouse.RData", sep=""))

bp.categories=GOdata[["biological_process"]][["categories"]]
bp.genelist=GOdata[["biological_process"]][["genelist"]]

################################################################################

network=read.table(paste(pathResults, type, "/common_interactions_FDR",maxFDR,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

################################################################################

tpm=read.table(paste(pathDatasets, "SupplementaryDataset6/NormalizedTPM_OrthoGenes_Mouse_Rat.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

rownames(tpm)=tpm$ID

info=tpm[,c(1,2)]
tpm=tpm[,-c(1,2)]

lnc=info$ID[which(info$GeneType=="lncRNA")]
pc=info$ID[which(info$GeneType=="protein_coding")]
pc.mouse=unlist(lapply(pc, function(x) unlist(strsplit(x, split="_"))[1]))
names(pc.mouse)=pc

meantpm.mouse=read.table(paste(pathDatasets, "SupplementaryDataset3/Statistics_AverageTPM_Mouse.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(meantpm.mouse)=meantpm.mouse$GeneID

meantpm.rat=read.table(paste(pathDatasets, "SupplementaryDataset3/Statistics_AverageTPM_Rat.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(meantpm.rat)=meantpm.rat$GeneID

################################################################################

results=data.frame("ID"=character(0), "NbPCConnections"=integer(0),"GOAccession"=character(0), "GOName"=character(0), "Enrichment"=numeric(0), "FDR"=numeric(0),  stringsAsFactors=F)

for(i in 1:length(lnc)){
  print(i)
  
  gene=lnc[i]

  w=which(network$ID1==gene | network$ID2==gene)
  
  if(length(w)>=minConnections){
    pcconn=intersect(c(network$ID1[w], network$ID2[w]), pc)
    
    if(length(pcconn)>=minConnections){
      pcconn=pc.mouse[pcconn]

      go=compute.go.enrichment(pcconn, pc.mouse, bp.categories, bp.genelist)

      if(length(which(go$FDR<maxGOFDR))>0){
        go=go[which(go$FDR<maxGOFDR),]
        go=go[order(go$FDR),]
        this.res=data.frame("ID"=rep(gene, dim(go)[1]), "NbPCConnections"=rep(length(pcconn), dim(go)[1]), "GOAccession"=go$GOAccession, "GOName"=go$GOName, "Enrichment"=go$Enrichment, "FDR"=go$FDR, stringsAsFactors=F)

        results=rbind(results, this.res)
        
      }
    }
  }
}

results$ID.Mouse=unlist(lapply(results$ID, function(x) unlist(strsplit(x, split="_"))[1]))
results$ID.Rat=unlist(lapply(results$ID, function(x) unlist(strsplit(x, split="_"))[2]))

results$MaxSample.Mouse=meantpm.mouse[results$ID.Mouse,"MaxSample"]

results$ExpressionSpecificity.Mouse=meantpm.mouse[results$ID.Mouse,"ExpressionSpecificity"]
results$MaxSample.Rat=meantpm.rat[results$ID.Rat,"MaxSample"]
results$ExpressionSpecificity.Rat=meantpm.rat[results$ID.Rat,"ExpressionSpecificity"]

results=results[,c("ID","ID.Mouse", "ID.Rat", setdiff(colnames(results), c("ID","ID.Mouse", "ID.Rat")))]

write.table(results, file=paste(pathResults, type, "/GOEnrichment_BiologicalProcess_lncRNAs_min", minConnections,"PCConnections.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

################################################################################



