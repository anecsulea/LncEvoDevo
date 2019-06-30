################################################################

path="LncEvoDevo/"
pathCodingPotential=paste(path, "results/coding_potential/", sep="")
pathAnnot=paste(path, "data/ensembl_annotations/", sep="")

#################################################################

release=94

prefixes=c(paste("FilteredTranscripts_StringTie_Ensembl",release,sep=""))
names(prefixes)=c("StringTie")

#################################################################

for(sp in c("Mouse", "Rat", "Chicken")){
  for(annot in c("StringTie")){
    class=read.table(paste(pathCodingPotential, sp, "/CombinedGeneClassification_", prefixes[annot],".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

    print(paste(sp, annot))

    pc.correct=class$GeneID[which(class$GeneBiotype=="protein_coding" & class$Class=="coding")]
    pc.wrong=class$GeneID[which(class$GeneBiotype=="protein_coding" & class$Class=="noncoding")]
    prop.correct.pc=length(pc.correct)/(length(pc.correct)+length(pc.wrong))
    
    linc.correct=class$GeneID[which(class$GeneBiotype=="lincRNA" & class$Class=="noncoding")]
    linc.wrong=class$GeneID[which(class$GeneBiotype=="lincRNA" & class$Class=="coding")]
    prop.correct.linc=length(linc.correct)/(length(linc.correct)+length(linc.wrong))
    
    mir.correct=class$GeneID[which(class$GeneBiotype=="miRNA" & class$Class=="noncoding")]
    mir.wrong=class$GeneID[which(class$GeneBiotype=="miRNA" & class$Class=="coding")]
    prop.correct.mir=length(mir.correct)/(length(mir.correct)+length(mir.wrong))


    snorna.correct=class$GeneID[which(class$GeneBiotype=="snoRNA" & class$Class=="noncoding")]
    snorna.wrong=class$GeneID[which(class$GeneBiotype=="snoRNA" & class$Class=="coding")]
    prop.correct.snorna=length(snorna.correct)/(length(snorna.correct)+length(snorna.wrong))
    
    rrna.correct=class$GeneID[which(class$GeneBiotype=="rRNA" & class$Class=="noncoding")]
    rrna.wrong=class$GeneID[which(class$GeneBiotype=="rRNA" & class$Class=="coding")]
    prop.correct.rrna=length(rrna.correct)/(length(rrna.correct)+length(rrna.wrong))


    print(paste(round(100*prop.correct.pc, digits=2), "% correct for protein-coding genes"))
    print(paste(round(100*prop.correct.linc, digits=2), "% correct for lincRNA genes"))
    print(paste(round(100*prop.correct.mir, digits=2), "% correct for miRNA genes"))
    print(paste(round(100*prop.correct.snorna, digits=2), "% correct for snoRNA genes"))
    print(paste(round(100*prop.correct.rrna, digits=2), "% correct for rRNA genes"))
    
  }
}

#################################################################





