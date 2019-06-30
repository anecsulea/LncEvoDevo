####################################################################

source("parameters.R")

pathDatasets=paste(path, "supplementary_datasets/", sep="")

splist=c("Mouse", "Rat", "Chicken")

####################################################################

allfiles=system(paste("ls ", pathDatasets, "SupplementaryDataset2/", sep=""),intern=T)

objects=c("rawtpm", "normtpm", "kcounts", "readcounts", "downsampled")
filenames=c("KallistoRawPM", "KallistoNormalizedTPM", "KallistoEstimatedCounts", "UniqueReadCounts", "UniqueReadCounts_Downsampled")

####################################################################

for(sp in splist){

  for(i in 1:length(objects)){
    o=objects[i]
    f=filenames[i]

    if(i==5){
      ff=grep(sp, grep(f, allfiles, value=T), value=T)
    } else{
      ff=paste(f,"_",sp,".txt",sep="")
    }
     
    print(paste(o, f, ff))
    data=read.table(paste(pathDatasets, "SupplementaryDataset2/",ff, sep=""), h=T, stringsAsFactors=F, sep="\t")
    rownames(data)=data$GeneID
    data=data[,-which(colnames(data)=="GeneID")]

    assign(o, data, envir=.GlobalEnv)
  }
  
  save(list=objects, file=paste("RData/data.expression.",sp,".RData", sep=""))
}

####################################################################
