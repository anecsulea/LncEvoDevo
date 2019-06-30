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
  
  
  stats=read.table(paste(pathMatrices,type,"/",sp,"/stats_matrices.txt",sep=""), h=T, stringsAsFactors=F)
  
  ok=which(stats$PropSynonymous.CDS>=0.6 & stats$PropNonsense.CDS<0.01)
  
  kept.informants=stats$Species2[ok]
  
  writeLines(kept.informants, con=paste(pathMatrices,type,"/",sp,"/selected_informant_species.txt",sep=""))
  
}

#######################################################################################################################
