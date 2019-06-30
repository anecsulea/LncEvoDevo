########################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathPromoters=paste(path, "results/promoter_analysis/", sep="")
pathRepeats=paste(path, "results/overlap_repeats/", sep="")

promsize="1kb"

########################################################################

for(ref in c("Mouse", "Rat")){
  
  tg=setdiff(c("Mouse", "Rat"), ref)
  
  
  ref.specific=read.table(paste(pathDatasets, "SupplementaryDataset8/Candidate_",ref,"Specific_LncRNAs_Min100Reads.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  
########################################################################
  
  ref.repeats=read.table(paste(pathRepeats, ref,"/OverlapRepeats_TEFamily_BothStrands_Promoters",promsize,"_FilteredTranscripts_StringTie_Ensembl94.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  colnames(ref.repeats)[1]="ID"
  ref.repeats$GeneID=unlist(lapply(ref.repeats$ID, function(x) unlist(strsplit(x, split="_"))[1]))
  
  print(length(unique(ref.repeats$GeneID)))
  print(length(ref.repeats$GeneID))
  
  tg.repeats=read.table(paste(pathRepeats, tg,"/OverlapRepeats_TEFamily_BothStrands_Promoters",promsize,"_From",ref,"_To",tg,"_ExonBlocks_FilteredTranscripts_StringTie_Ensembl94_FilteredProjectedExons_Step2.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  colnames(tg.repeats)[1]="ID"
  tg.repeats$GeneID=unlist(lapply(tg.repeats$ID, function(x) unlist(strsplit(x, split="_"))[1]))
  
  print(length(unique(tg.repeats$GeneID)))
  print(length(tg.repeats$GeneID))


  common=intersect(tg.repeats$GeneID, ref.repeats$GeneID)
  tg.repeats=tg.repeats[which(tg.repeats$GeneID%in%common),]
  ref.repeats=ref.repeats[which(ref.repeats$GeneID%in%common),]
  
########################################################################
  
  ref.specific.repeats=ref.repeats[which(ref.repeats$GeneID%in%ref.specific$GeneID),]
  tg.specific.repeats=tg.repeats[which(tg.repeats$GeneID%in%ref.specific$GeneID),]
  
  ########################################################################
  
  selected.families=c("Alu", "B2", "B4", "CR1", "ERV1", "ERVK", "ERVL", "ERVL.MaLR", "L1", "L2", "MIR", "Simple_repeat", "Low_complexity")
  erv.families=c( "ERV1", "ERVK", "ERVL", "ERVL.MaLR")
  line.families=c("L1", "L2")
  
  prop.genes.ref=list()
  prop.genes.tg=list()

  prop.genes.ref.specific=list()
  prop.genes.tg.specific=list()
  
  for(fam in selected.families){
    prop.genes.ref[[fam]]=length(which(ref.repeats[,fam]>0))/dim(ref.repeats)[1]
    prop.genes.tg[[fam]]=length(which(tg.repeats[,fam]>0))/dim(tg.repeats)[1]

    prop.genes.ref.specific[[fam]]=length(which(ref.specific.repeats[,fam]>0))/dim(ref.specific.repeats)[1]
    prop.genes.tg.specific[[fam]]=length(which(tg.specific.repeats[,fam]>0))/dim(tg.specific.repeats)[1]
  }

  
  print(unlist(prop.genes.ref.specific)/unlist(prop.genes.tg.specific))
  print(unlist(prop.genes.ref)/unlist(prop.genes.tg))


########################################################################
}
