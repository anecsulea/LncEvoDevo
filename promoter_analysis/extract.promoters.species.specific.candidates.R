################################################################################

library(seqinr)

################################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathPromoters=paste(path, "results/promoter_analysis/", sep="")

release=94

################################################################################

splist=c("Mouse", "Rat")

for(minreads in c(100, 250)){
  for(ref in splist){

    tg=setdiff(splist, ref)

    candidates=read.table(paste(pathDatasets, "SupplementaryDataset8/Candidate_",ref,"Specific_LncRNAs_Min", minreads,"Reads.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    candidates$TSS=candidates$Start
    candidates$TSS[which(candidates$Strand==-1)]=candidates$End[which(candidates$Strand==-1)]
    candidates$PromoterStart=candidates$TSS
    candidates$PromoterStart[which(candidates$Strand==1)]=candidates$TSS[which(candidates$Strand==1)]-1-1000

    rownames(candidates)=candidates$GeneID
       
    coords.ref=read.table(paste(pathPromoters, ref,"/PromoterCoords_1kb_FilteredTranscripts_StringTie_Ensembl",release,".bed", sep=""), h=F, stringsAsFactors=F)
    colnames(coords.ref)=c("Chr", "Start", "End", "ID","Score", "Strand")
    coords.ref$GeneID=unlist(lapply(coords.ref[,4], function(x) unlist(strsplit(x, split="_"))[1]))
    coords.ref=coords.ref[which(coords.ref$GeneID%in%candidates$GeneID),]

    coords.ref=coords.ref[which(coords.ref$Start==candidates[coords.ref$GeneID,"PromoterStart"]),]

    print(all(candidates$GeneID%in%coords.ref$GeneID))


    print("reading reference sequences")
    fasta.ref=read.fasta(paste(pathPromoters,ref, "/PromoterSequences_1kb_FilteredTranscripts_StringTie_Ensembl",release,".fa", sep=""),seqtype="DNA", forceDNAtolower=FALSE)
    selected.ref=fasta.ref[which(names(fasta.ref)%in%coords.ref$ID)]
    write.fasta(selected.ref,names=names(selected.ref), file.out=paste(pathPromoters,ref, "/PromoterSequences_Candidate_",ref,"Specific_LncRNAs_Min", minreads,"Reads.fa",sep=""))
    

    print("reading tg sequences")
    fasta.tg=read.fasta(paste(pathPromoters, tg, "/From",ref,"_PromoterCoords_1kb_FilteredTranscripts_StringTie_Ensembl",release,".fa", sep=""), seqtype="DNA", forceDNAtolower=FALSE)
    selected.tg=fasta.tg[which(names(fasta.tg)%in%paste(tg, coords.ref$ID, sep="_"))]

    write.fasta(selected.tg,names=names(selected.tg), file.out=paste(pathPromoters,ref, "/PromoterSequences_Candidate_",ref,"Specific_LncRNAs_Min", minreads,"Reads.fa",sep=""))

    stop()
    
  }
}

################################################################################



