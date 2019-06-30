###############################################################################

path="LncEvoDevo/"
pathCodingPotential=paste(path, "results/coding_potential/", sep="")

annot="FilteredTranscripts_StringTie_Ensembl89"

aln.list=list()
aln.list[["Mouse"]]=c("Mouse_60way")
aln.list[["Rat"]]=c("Human_100way", "Rat_20way")
aln.list[["Chicken"]]=c("Human_100way")


###############################################################################

for(sp in c("Mouse", "Rat", "Chicken")){

  print(sp)

  for(aln in aln.list[[sp]]){
    print(aln)
    
    class.allexons=read.table(paste(pathCodingPotential, sp, "/GeneClassification_AllExons_CSFScores_",aln,"_windowSize75_", annot, ".txt", sep=""), h=T, stringsAsFactors=F)
    rownames(class.allexons)=class.allexons$GeneID
    
    class.noov=read.table(paste(pathCodingPotential, sp, "/GeneClassification_NonOverlappingRegions_CSFScores_",aln,"_windowSize75_", annot, ".txt", sep=""), h=T, stringsAsFactors=F)
    rownames(class.noov)=class.noov$GeneID
    
    only.allexons=setdiff(rownames(class.allexons), rownames(class.noov))
    only.noov=setdiff(rownames(class.noov), rownames(class.allexons))
    common.genes=intersect(rownames(class.allexons), rownames(class.noov))
    class.allexons.common=class.allexons[common.genes,]
    class.noov.common=class.noov[common.genes,]
    
    coding.noov=which(class.allexons.common$Class!=class.noov.common$Class & class.noov.common$Class=="coding")
    noncoding.noov=which(class.allexons.common$Class!=class.noov.common$Class & class.noov.common$Class=="noncoding")
    
    nb.agreement=length(which(class.allexons.common$Class==class.noov.common$Class))
    nb.disagreement=length(which(class.allexons.common$Class!=class.noov.common$Class))
    pc.agreement=round(100*nb.agreement/(nb.agreement+nb.disagreement), digits=2)
    
    print(paste(length(only.allexons), "genes present only in the all exon classification"))
    print(paste(length(only.noov), "genes present only in the non-overlapping exon classification"))
    
    print(paste(nb.agreement, " genes with the same classification (",pc.agreement,"%)", sep=""))
    print(paste(nb.disagreement, "genes with different classifications"))
    print(paste(length(coding.noov), "genes coding only with non-overlapping exons"))
    print(paste(length(noncoding.noov), "genes noncoding only with non-overlapping exons"))
  }
}

###############################################################################


