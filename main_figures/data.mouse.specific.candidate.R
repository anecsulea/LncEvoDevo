##########################################################

source("parameters.R")
pathGenomeBrowser=paste(path, "results/genome_browser/", sep="")
pathStringTie=paste(path, "results/stringtie_assembly/", sep="")
pathEnhancers=paste(path, "data/enhancers/Encode_YueLab/", sep="")

##########################################################

mouse.id="ENSMUSG00000109311"
mouse.chr="7"
mouse.start=89426306 -30000
mouse.end=89441928 +5000
mouse.strand="forward"
proj.mouse.start=89417975
proj.mouse.end=89423479

rat.id=NA
rat.chr="1"
proj.rat.start=153610811
proj.rat.end=153626049
rat.start=proj.rat.start-30000
rat.end=proj.rat.end+5000
rat.strand="forward"

##########################################################

release=94

annot.mouse=read.table(paste(pathStringTie, "Mouse/combined/ExonBlocks_FilteredTranscripts_StringTie_Ensembl",release,".txt",sep=""), h=F, stringsAsFactors=F)
annot.rat=read.table(paste(pathStringTie, "Rat/combined/ExonBlocks_FilteredTranscripts_StringTie_Ensembl",release,".txt",sep=""), h=F, stringsAsFactors=F)  

##########################################################

mouse.enhancers=read.table(paste(pathEnhancers,"predicted_enhancer_mouse/all_enhancers_500bp_mm10.bed",sep=""), h=F, stringsAsFactors=F)
mouse.enhancers=mouse.enhancers[which(mouse.enhancers[,1]==paste("chr",mouse.chr,sep="")),]

##########################################################

cov.mouse=list()
cov.rat=list()

for(sample in c("Kidney_Adult", "Kidney_Aged")){

  this.cov.mouse=read.table(paste(pathGenomeBrowser, "Mouse/", sample,"/coverage_unique_forward_",mouse.chr,".bedGraph", sep=""), h=F,stringsAsFactors=F)
  this.cov.mouse=this.cov.mouse[which(this.cov.mouse[,2]>=mouse.start & this.cov.mouse[,3]<=mouse.end),]

  cov.mouse[[sample]]=this.cov.mouse

  this.cov.rat=read.table(paste(pathGenomeBrowser, "Rat/", sample,"/coverage_unique_forward_", rat.chr,".bedGraph", sep=""), h=F,stringsAsFactors=F)
  this.cov.rat=this.cov.rat[which(this.cov.rat[,2]>=rat.start & this.cov.rat[,3]<=rat.end),]

  cov.rat[[sample]]=this.cov.rat
}

save(list=c("mouse.id","mouse.chr", "mouse.start", "mouse.end", "mouse.strand", "rat.chr", "rat.start", "rat.end", "rat.strand", "mouse.enhancers", "annot.mouse", "annot.rat", "cov.mouse", "cov.rat", "proj.rat.start", "proj.rat.end", "proj.mouse.start", "proj.mouse.end"), file="RData/data.mouse.specific.candidate.RData")

##########################################################
