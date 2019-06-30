#################################################################

path="LncEvoDevo/"
pathEnsembl=paste(path, "data/ensembl_ortho/", sep="")
pathAnnot=paste(path, "data/ensembl_annotations/", sep="")
pathWGA=paste(path, "results/ortho_genes/whole_genome_alignments/", sep="")

#################################################################

splist=list(c("Mouse", "Rat"), c("Mouse", "Chicken"), c("Rat", "Chicken"))
release=94

#################################################################

for(i in 1:3){
  
  sp1=splist[[i]][1]
  sp2=splist[[i]][2]
  
  info.sp1=read.table(paste(pathAnnot, sp1, "/GeneInfo_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  info.sp2=read.table(paste(pathAnnot, sp2, "/GeneInfo_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    
  lnc.sp1=info.sp1[which(info.sp1$biotype%in%c("lincRNA","antisense", "processed_transcript")),1]
  lnc.sp2=info.sp2[which(info.sp2$biotype%in%c("lincRNA","antisense", "processed_transcript")),1]
  
#################################################################
  
  ens=read.table(paste(pathEnsembl, "GeneFamilies_1to1_",sp1, sp2,"_ProteinCoding_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F)
  rownames(ens)=ens[,sp1]
  
  wga=read.table(paste(pathWGA, "ReciprocalBestHits_",sp1,"_",sp2,"_StringTie.txt",sep=""), h=T, stringsAsFactors=F)
  rownames(wga)=wga[,paste("ID",sp1,sep=".")]
  
#################################################################
  
  all.sp1=unique(c(ens[,sp1], wga[,paste("ID",sp1, sep=".")]))
  
  all.ortho=data.frame("IDRef"=all.sp1, "EnsemblOrtho"=rep(NA, length(all.sp1)), "WGAOrtho"=rep(NA, length(all.sp1)), stringsAsFactors=F)
  all.ortho$EnsemblOrtho[which(all.ortho$IDRef%in%ens[,sp1])]=ens[all.ortho[which(all.ortho$IDRef%in%ens[,sp1]), "IDRef"],sp2]

  all.ortho$WGAOrtho[which(all.ortho[,"IDRef"]%in%wga[,paste("ID",sp1, sep=".")])]=wga[all.ortho[which(all.ortho$IDRef%in%wga[,paste("ID",sp1, sep=".")]), "IDRef"],paste("ID",sp2, sep=".")]
                      
  #################################################################
  
  notfound=all.ortho[which(is.na(all.ortho$WGAOrtho)),]
  different=all.ortho[which(all.ortho$EnsemblOrtho!=all.ortho$WGAOrtho),]
  same=all.ortho[which(all.ortho$EnsemblOrtho==all.ortho$WGAOrtho),]
  lnc=all.ortho[which(all.ortho$IDRef%in%lnc.sp1 | all.ortho$WGAOrtho%in%lnc.sp2),]

  nb.notfound=dim(notfound)[1]
  nb.different=dim(different)[1]
  nb.same=dim(same)[1]
  nb.lnc=dim(lnc)[1]

  pcok=round(100*nb.same/dim(ens)[1], digits=1)

  print(paste(sp1, sp2))

  print(paste(nb.notfound, "ortho not found"))
  print(paste(nb.different, "ortho different"))
  print(paste(nb.same, "ortho identical"))
  print(paste(pcok, "% correct"))
  print(paste(nb.lnc, "ortho lncRNA"))
 
  print("")
}

#################################################################
