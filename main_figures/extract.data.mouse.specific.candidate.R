##########################################################

source("parameters.R")
pathGenomeBrowser=paste(path, "results/genome_browser/", sep="")

##########################################################

mouse.id="ENSMUSG00000109311"
mouse.chr="7"
mouse.start=89426306 -20000
mouse.end=89441928 +20000
mouse.strand="forward"


rat.id=NA
rat.chr="1"
rat.start=153610811-20000
rat.end=153626049+10000
rat.strand="forward"

samples=c("Kidney_Adult", "Kidney_Aged")

for(sample in samples){
  sp="Mouse"
  this.chr=mouse.chr
  this.strand=mouse.strand
  
  system(paste("zcat ", pathGenomeBrowser, sp, "/", sample, "/coverage_unique_", this.strand,".bedGraph.gz | grep ^",this.chr," > ",pathGenomeBrowser, sp, "/", sample, "/coverage_unique_",this.strand,"_",this.chr,".bedGraph",sep=""))

  sp="Rat"
  this.chr=rat.chr
  this.strand=rat.strand
  
  system(paste("zcat ", pathGenomeBrowser, sp, "/", sample, "/coverage_unique_", this.strand,".bedGraph.gz | grep ^",this.chr," > ",pathGenomeBrowser, sp, "/", sample, "/coverage_unique_",this.strand,"_",this.chr,".bedGraph",sep=""))
  
}

##########################################################
