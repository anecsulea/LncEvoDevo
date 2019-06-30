########################################################################

path="LncEvoDevo/"
pathCSF=paste(path, "results/coding_potential/",sep="")
pathEnsembl=paste(path, "data/ensembl_annotations/", sep="")

########################################################################

aln=list()
aln[["Mouse"]]=c("Mouse_60way")
aln[["Chicken"]]=c("Human_100way")
aln[["Rat"]]=c("Human_100way", "Rat_20way")

########################################################################

for(sp in c("Mouse", "Rat", "Chicken")){
  for(this.aln in aln[[sp]]){
    geneinfo=read.table(paste(pathEnsembl, sp, "/GeneInfo_Ensembl89.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    pcgenes=geneinfo$stable_id[which(geneinfo$biotype=="protein_coding" & geneinfo$status=="KNOWN")]
    
    cds=read.table(paste(pathEnsembl, sp, "/CDSCoords_Ensembl89.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    cds=cds[which((!is.na(cds$Genomic.coding.start)) & (!is.na(cds$Genomic.coding.end))),]
    cds$ID=paste(cds$Chromosome.scaffold.name, cds$Genomic.coding.start,  cds$Genomic.coding.end,  cds$Strand, sep=",")
    cds=cds[which(cds$Gene.stable.ID%in%pcgenes),]
    
    csf=read.table(paste(pathCSF, sp, "/OverlapExons_CSFScores_",this.aln,"_windowSize75_FilteredTranscripts_Ensembl89.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    csf$ID=paste(csf$Chr, csf$Start, csf$End, csf$Strand, sep=",")
    csf=csf[which(csf$TotalLength==csf$LengthCoveredCSF),]
    
    common=intersect(cds$ID, csf$ID)
    
    csf.common=csf[which(csf$ID%in%common),]
    
    max.csf.gene=tapply(csf.common$MaxScore, as.factor(csf.common$GeneID), max)
    
    stop()
  }
}

########################################################################
