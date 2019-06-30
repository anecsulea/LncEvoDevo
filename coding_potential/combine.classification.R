#################################################################

path="LncEvoDevo/"
pathCodingPotential=paste(path, "results/coding_potential/", sep="")
pathAnnot=paste(path, "data/ensembl_annotations/", sep="")

#################################################################

release=94

csfaln=list()
csfaln[["Mouse"]]=c("Mouse_60way")
csfaln[["Rat"]]=c("Rat_20way", "Human_100way")
csfaln[["Chicken"]]=c("Human_100way")

prefixes=c(paste("FilteredTranscripts_StringTie_Ensembl",release, sep=""))
names(prefixes)=c("StringTie")

#################################################################

for(sp in c("Mouse", "Rat", "Chicken")){
  print(sp)
  
  geneinfo=read.table(paste(pathAnnot,sp, "/GeneInfo_Ensembl",release,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(geneinfo)=geneinfo$stable_id
  
  for(annot in c("StringTie")){
    print(annot)

    ## BlastX
    
    blastx.swiss=read.table(paste(pathCodingPotential, sp, "/GeneClassification_", prefixes[annot], "_rm_vs_SwissProt.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    rownames(blastx.swiss)=blastx.swiss$GeneID

    blastx.pfam=read.table(paste(pathCodingPotential, sp, "/GeneClassification_", prefixes[annot], "_rm_vs_Pfam.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    rownames(blastx.pfam)=blastx.pfam$GeneID

    ## CSF

    all.class=data.frame("GeneID"=blastx.swiss$GeneID, "BlastX.Swiss"=blastx.swiss$Class, "BlastX.Pfam"=blastx.pfam[blastx.swiss$GeneID,"Class"], stringsAsFactors=F)
    
    for(aln in csfaln[[sp]]){
      print("using non-overlapping regions for CSF")
      
      this.csf=read.table(paste(pathCodingPotential, sp, "/GeneClassification_NonOverlappingRegions_CSFScores_", aln,"_windowSize75_", prefixes[annot],".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

      rownames(this.csf)=this.csf$GeneID
      
      all.class[,paste("CSF.",aln,sep="")]=this.csf[all.class$GeneID,"Class"]
    }


    
    ## final classification
  
    final.class=rep("NA", dim(all.class)[1])
    final.class[which(apply(all.class,1,function(x) length(which(x=="coding")))>0)]="coding"
    final.class[which(apply(all.class,1,function(x) length(which(x=="coding")))==0 & apply(all.class,1,function(x) length(which(x=="noncoding")))>0)]="noncoding"

    all.class$Class=final.class

    all.class$GeneBiotype=rep(NA, dim(all.class)[1])
    all.class[which(all.class$GeneID%in%rownames(geneinfo)),"GeneBiotype"]=geneinfo[all.class[which(all.class$GeneID%in%rownames(geneinfo)),"GeneID"],"biotype"]
    
    write.table(all.class, file=paste(pathCodingPotential, sp, "/CombinedGeneClassification_", prefixes[annot],".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  }
}

#################################################################
