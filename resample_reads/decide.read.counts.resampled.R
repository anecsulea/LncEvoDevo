#######################################################################

path="LncEvoDevo/"
pathHisat=paste(path, "results/hisat/", sep="")

options(scipen=999)

otherstrain=list()
otherstrain[["Mouse"]]=c(NA)
otherstrain[["Rat"]]=c("Brain_Adult3", "Brain_Adult4", "Kidney_Adult3", "Kidney_Adult4", "Liver_Adult3", "Liver_Adult4", "Testis_Adult4")
otherstrain[["Chicken"]]=c(NA)

nbresampled=list()
nbresampled[["Mouse"]]=6e7
nbresampled[["Rat"]]=6e7
nbresampled[["Chicken"]]=5e7

accepted.tissage=list()
accepted.tissage[["Chicken"]]=c("Brain_EarlyEmbryo", "Brain_LateEmbryo", "Kidney_EarlyEmbryo", "Kidney_LateEmbryo", "Liver_EarlyEmbryo", "Liver_LateEmbryo")

#######################################################################

for(sp in c("Mouse", "Rat", "Chicken")){
  print(sp)
  
  table.nbreads=read.table(paste(pathHisat, sp, "/nb_unique_reads_noMT.txt", sep=""), h=F, stringsAsFactors=F)
  table.nbreads=table.nbreads[which(!table.nbreads[,1]%in%otherstrain[[sp]]),]

  nbreads=table.nbreads[,2]
  samples=table.nbreads[,1]
  
  names(nbreads)=samples
  tissues=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))
  ages=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[2]))
  ages=unlist(lapply(ages, function(x) substr(x, 1, nchar(x)-1)))

  tissage=paste(tissues, ages, sep="_")
  nbrepli=as.numeric(table(tissage))
  names(nbrepli)=levels(as.factor(tissage))

  totreads=tapply(nbreads, as.factor(tissage), sum)
  names(totreads)=levels(as.factor(tissage))

  all.read.stats=data.frame("Sample"=samples, "TissueAge"=tissage, "NbReads"=nbreads, stringsAsFactors=F)
 
  combined.read.stats=data.frame("TissueAge"=levels(as.factor(tissage)), "TotReads"=totreads, "NbReplicates"=nbrepli, stringsAsFactors=F)
  rownames(combined.read.stats)=combined.read.stats$TissueAge

  all.read.stats$NbResampled=nbresampled[[sp]]/combined.read.stats[all.read.stats$TissueAge, "NbReplicates"]

  all.read.stats$TissueAge[which(all.read.stats$Sample%in%c("Sertoli", "Spermatogonia", "Spermatocytes", "Spermatids", "Spermatozoa"))]=all.read.stats$Sample[which(all.read.stats$Sample%in%c("Sertoli", "Spermatogonia", "Spermatocytes", "Spermatids", "Spermatozoa"))]

  if(sp=="Chicken"){
    all.read.stats=all.read.stats[which(all.read.stats$TissueAge%in%accepted.tissage[[sp]]),]
  }
  
  if(sp=="Mouse"){
    all.read.stats[which(all.read.stats$Sample=="Liver_Adult2"), "NbResampled"]=25e6
    all.read.stats[which(all.read.stats$Sample=="Liver_Adult1"), "NbResampled"]=35e6
  }

  if(sp=="Rat"){
    all.read.stats[which(all.read.stats$Sample=="Kidney_Aged1"), "NbResampled"]=32e6
    all.read.stats[which(all.read.stats$Sample=="Kidney_Aged2"), "NbResampled"]=28e6
  }

  if(any(all.read.stats$NbResampled>all.read.stats$NbReads)){
    print("weird! saw resampled read counts larger than total read counts")
  }
  
  write.table(all.read.stats, file=paste(pathHisat, sp, "/nb_resampled_unique_reads_noMT.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

}

#######################################################################


