#########################################################################

load("RData/data.expression.ortho.RData") 

source("parameters.R")

#########################################################################

compute.euclidean.distance <- function(expdata, minexp, discarded.samples=NA, norm.method="max"){
  if(any(!is.na(discarded.samples))){
    notok=unlist(lapply(discarded.samples, function(x) grep(x, colnames(expdata), value=T)))
    notok=unique(notok)

    if(length(notok)>0){
      print("discarding following samples:")
      print(notok)
      
      expdata=expdata[,-which(colnames(expdata)%in%notok)]
    }
  }
  
  this.samples=colnames(expdata)
  this.species=unlist(lapply(this.samples, function(x) unlist(strsplit(x, split="_"))[1]))
  this.tissues=unlist(lapply(this.samples, function(x) unlist(strsplit(x, split="_"))[2]))
  this.ages=unlist(lapply(this.samples, function(x) unlist(strsplit(x, split="_"))[3]))
  this.ages=unlist(lapply(this.ages, function(x) substr(x, 1, nchar(x)-1)))

  this.tissage=paste(this.tissues, this.ages, sep="_")
  this.sptissage=as.factor(paste(this.species, this.tissues, this.ages, sep="_"))

  this.meanexp=t(apply(expdata, 1, function(x) tapply(as.numeric(x), this.sptissage, mean)))
  
  sp1=unique(this.species)[1]
  sp2=unique(this.species)[2]
  
  expdata1=expdata[,which(this.species==sp1)]
  expdata2=expdata[,which(this.species==sp2)]

  this.tissage1=as.factor(this.tissage[which(this.species==sp1)])
  this.tissage2=as.factor(this.tissage[which(this.species==sp2)])

  mean.expdata1=t(apply(expdata1, 1, function(x) tapply(as.numeric(x), this.tissage1, mean)))
  mean.expdata2=t(apply(expdata2, 1, function(x) tapply(as.numeric(x), this.tissage2, mean)))

  rownames(mean.expdata1)=rownames(expdata)
  rownames(mean.expdata2)=rownames(expdata)
  
  ok=rownames(expdata)[which(apply(mean.expdata1, 1, max)>minexp | apply(mean.expdata2, 1, max)>minexp)]

  print(paste(length(ok), "genes with max expression greater than",minexp))
  
  mean.expdata1=mean.expdata1[ok,]
  mean.expdata2=mean.expdata2[ok,]
  
  nbgenes=length(ok)

   if(norm.method=="none"){
    print("no normalization, computing distance on raw data")
  } else{
    if(norm.method=="max"){
      print("dividing by maximum expression")
      mean.expdata1=mean.expdata1/apply(mean.expdata1,1,max)
      mean.expdata2=mean.expdata2/apply(mean.expdata2,1,max)
    } else{
      if(norm.method=="sum"){
        print("dividing by sum expression")
        mean.expdata1=mean.expdata1/apply(mean.expdata1,1,sum)
        mean.expdata2=mean.expdata2/apply(mean.expdata2,1,sum)
      } else{
        if(norm.method=="z-score" | norm.method=="zscore" ){
          print("computing z-score")
          mean.expdata1=t(apply(mean.expdata1, 1, function(x) (x-mean(x))/sd(x)))
          mean.expdata2=t(apply(mean.expdata2, 1, function(x) (x-mean(x))/sd(x)))
        } else{
          stop("unkown normalization method")
        }
      }
    }
  }
  
  distance=unlist(lapply(1:nbgenes, function(x) dist(matrix(c(as.numeric(mean.expdata1[x,]), as.numeric(mean.expdata2[x,])), nrow=2, byrow=T))[1]))
  names(distance)=ok

  return(list("distance"=distance, "norm.exp.1"=mean.expdata1, "norm.exp.2"=mean.expdata2))
}

########################################################################
########################################################################

minexp=0

tpm.orthopc=exportho.mr[which(exportho.mr$GeneType=="protein_coding"), which(!colnames(exportho.mr)%in%c("ID", "GeneType"))]
tpm.ortholnc=exportho.mr[which(exportho.mr$GeneType=="lncRNA"), which(!colnames(exportho.mr)%in%c("ID", "GeneType"))]

rownames(tpm.orthopc)=exportho.mr$ID[which(exportho.mr$GeneType=="protein_coding")]
rownames(tpm.ortholnc)=exportho.mr$ID[which(exportho.mr$GeneType=="lncRNA")]

hasnapc=apply(tpm.orthopc, 1, function(x) any(is.na(x)))
tpm.orthopc=tpm.orthopc[which(!hasnapc),]

hasnalnc=apply(tpm.ortholnc, 1, function(x) any(is.na(x)))
tpm.ortholnc=tpm.ortholnc[which(!hasnalnc),]

exptype="tpm"

########################################################################

for(age in age.order){
  for(genetype in c("pc", "lnc")){
    expdata=get(paste(exptype, ".ortho", genetype, sep=""))
    expdata=expdata[,grep(age, colnames(expdata))]
    
    for(method in c("sum")){
      results.allsamples=compute.euclidean.distance(expdata, minexp, norm.method=method)
      
      dist.allsamples=results.allsamples[["distance"]]
      
      expdiv=data.frame("GeneID"=names(dist.allsamples), "ED.AllSamples"=dist.allsamples, stringsAsFactors=F)
      
      save(expdiv, file=paste("RData/data.expression.divergence.",exptype,".", genetype,".",method,".",age,".RData", sep=""))
    }
  }
}

########################################################################

