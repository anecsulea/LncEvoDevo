#######################################################################################################################

path="LncEvoDevo/"
pathMatrices=paste(path, "results/CSF/matrices/", sep="")

#######################################################################################################################

species=c("Mouse", "Rat", "Chicken", "Human", "Rat")
types=c("Mouse_60way", "Human_100way", "Human_100way", "Human_100way", "Rat_20way")
refs=c("mm10", "rn6", "galGal4", "hg38", "rn6")

#######################################################################################################################

for(i in c(5)){ ## 1:5
  sp=species[i]
  type=types[i]
  ref=refs[i]
  
  pathCDS=paste(pathMatrices,type,"/",sp,"/matrices_CDS.txt",sep="")
  pathIntrons=paste(pathMatrices,type,"/",sp,"/matrices_introns.txt",sep="")
  
  allsp=system(paste("grep ^# ",pathCDS,sep=""), intern=T)
  tgsp=unlist(lapply(allsp, function(x) unlist(strsplit(x,split=" "))[3]))
  
#######################################################################################################################
  
  nbcodons=64*63
  
  ## output
  
  outfile=file(paste(pathMatrices,type,"/",sp,"/stats_matrices.txt",sep=""), open="w")
  
  writeLines("Species1\tSpecies2\tNbSynonymous.CDS\tNbMissense.CDS\tNbNonsense.CDS\tPropSynonymous.CDS\tPropNonsense.CDS\tNbSynonymous.Introns\tNbMissense.Introns\tNbNonsense.Introns\tPropSynonymous.Introns\tPropNonsense.Introns", outfile)
  
  for(tg in tgsp){
    
    print(paste(ref,tg))
    
    system(paste("grep -A ",nbcodons," \"# ",ref," ",tg,"\" ",pathCDS," > tmp_CDS_",ref,"_",tg,".txt",sep=""))
    system(paste("grep -A ",nbcodons," \"# ",ref," ",tg,"\" ",pathIntrons," > tmp_Introns_",ref,"_",tg,".txt",sep=""))
    
    tmp=read.table(paste("tmp_CDS_",ref,"_",tg,".txt",sep=""), stringsAsFactors=F)
    
    rownames(tmp)=paste(tmp$V1,tmp$V3,sep="_")
    
    syn.cds=sum(tmp$V5[which(tmp$V2==tmp$V4 & tmp$V2!="*" & tmp$V4!="*")])
    missense.cds=sum(tmp$V5[which(tmp$V2!=tmp$V4 & tmp$V2!="*" & tmp$V4!="*")])
    nonsense.cds=sum(tmp$V5[which(tmp$V2!=tmp$V4 & (tmp$V2=="*" | tmp$V4=="*"))])
    
    propsyn.cds=round(syn.cds/(syn.cds+missense.cds+nonsense.cds),digits=2)
    propnon.cds=round(nonsense.cds/(syn.cds+missense.cds+nonsense.cds),digits=2)
    
    system(paste("rm tmp_CDS_",ref,"_",tg,".txt",sep=""))
    
    tmp=read.table(paste("tmp_Introns_",ref,"_",tg,".txt",sep=""), stringsAsFactors=F)
    
    rownames(tmp)=paste(tmp$V1,tmp$V3,sep="_")
    
    syn.int=sum(tmp$V5[which(tmp$V2==tmp$V4 & tmp$V2!="*" & tmp$V4!="*")])
    missense.int=sum(tmp$V5[which(tmp$V2!=tmp$V4 & tmp$V2!="*" & tmp$V4!="*")])
    nonsense.int=sum(tmp$V5[which(tmp$V2!=tmp$V4 & (tmp$V2=="*" | tmp$V4=="*"))])
    
    propsyn.int=round(syn.int/(syn.int+missense.int+nonsense.int),digits=2)
    propnon.int=round(nonsense.int/(syn.int+missense.int+nonsense.int),digits=2)
    
    system(paste("rm tmp_Introns_",ref,"_",tg,".txt",sep=""))
    
    writeLines(paste(ref, tg, syn.cds, missense.cds, nonsense.cds, propsyn.cds, propnon.cds, syn.int, missense.int, nonsense.int, propsyn.int, propnon.int, sep="\t"), outfile) 
    
  }
  
  close(outfile)

}

#######################################################################################################################
