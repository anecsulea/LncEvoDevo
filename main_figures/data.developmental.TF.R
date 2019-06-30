##########################################################################

load("RData/data.gene.ontology.Mouse.RData")
load("RData/data.annotations.Mouse.RData")

##########################################################################

bp=GOdata[["biological_process"]]

bp.genelist=bp[["genelist"]]
bp.categories=bp[["categories"]]
  
##########################################################################

dev=bp.categories$GO.term.accession[which(bp.categories$GO.term.name%in%c("multicellular organism development", "system development","embryonic organ development", "animal organ development", "pattern specification process"))]

exp=bp.categories$GO.term.accession[which(bp.categories$GO.term.name%in%c("regulation of gene expression", "gene expression", "positive regulation of gene expression", "negative regulation of gene expression", "regulation of transcription by RNA polymerase II", "regulation of transcription, DNA-templated"))]


dev.tf=intersect(unlist(bp.genelist[dev]), unlist(bp.genelist[exp]))
  
##########################################################################

save(dev.tf, file="RData/data.developmental.TF.RData")

##########################################################################
