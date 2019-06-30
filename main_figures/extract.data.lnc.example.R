############################################################

if(load){
  load("RData/data.annotations.Chicken.RData")
  all.info.chicken=allinfo
  
  load("RData/data.annotations.Rat.RData")
  all.info.rat=allinfo
  
  load("RData/data.annotations.Mouse.RData")
  all.info.mouse=allinfo

  source("parameters.R")
  
  pathGenomeBrowser=paste(path, "results/genome_browser/", sep="")
  pathAnnot=paste(path, "results/stringtie_assembly/", sep="")

  load=F
}

############################################################

## #genes.mouse=c("Locus_31563", "ENSMUSG00000021767"),
## genes.mouse=c("ENSMUSG00000087211", "ENSMUSG00000018698")
## #genes.rat=c("Locus_44367", "ENSRNOG00000012850")
## genes.rat=c("ENSRNOG00000055694", "ENSRNOG00000002812")
## #genes.chicken=c("Locus_15781", "ENSGALG00000005035")
## genes.chicken=c("Locus_1671", "ENSGALG00000005409")

genes.mouse=c("MSTRG.15799")
genes.rat=c("MSTRG.66441")
genes.chicken=c("MSTRG.32907")

chr.mouse=all.info.mouse[genes.mouse[1],"Chr"]
chr.rat=all.info.rat[genes.rat[1],"Chr"]
chr.chicken=all.info.chicken[genes.chicken[1],"Chr"]

############################################################

for(sp in c("Mouse", "Rat", "Chicken")){
  this.chr=get(paste("chr", tolower(sp), sep="."))
  
  for(sample in c("Brain_EarlyEmbryo", "Brain_LateEmbryo", "Kidney_EarlyEmbryo", "Kidney_LateEmbryo")){
    print(paste(sp, sample))
    
    system(paste("zcat ", pathGenomeBrowser, sp, "/", sample, "/coverage_unique_forward.bedGraph.gz | grep ^",this.chr," > ",pathGenomeBrowser, sp, "/", sample, "/coverage_unique_forward_",this.chr,".bedGraph",sep=""))
    
    system(paste("zcat ", pathGenomeBrowser, sp, "/", sample, "/coverage_unique_reverse.bedGraph.gz | grep ^",this.chr," > ",pathGenomeBrowser, sp, "/", sample, "/coverage_unique_reverse_",this.chr,".bedGraph",sep=""))
  }
}

############################################################

for(gene in genes.mouse){
  system(paste("grep ",gene," ",pathAnnot, "Mouse/combined/ExonBlocks_FilteredTranscripts_StringTie_Ensembl94.txt > lnc_examples/Mouse_", gene, ".txt",sep=""))
}

for(gene in genes.rat){
  system(paste("grep ",gene," ",pathAnnot, "Rat/combined/ExonBlocks_FilteredTranscripts_StringTie_Ensembl94.txt > lnc_examples/Rat_", gene, ".txt",sep=""))
}

for(gene in genes.chicken){
  system(paste("grep ",gene," ",pathAnnot, "Chicken/combined/ExonBlocks_FilteredTranscripts_StringTie_Ensembl94.txt > lnc_examples/Chicken_", gene, ".txt",sep=""))
}

############################################################

