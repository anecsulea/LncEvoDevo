########################################################################

set.seed(19)


minTPM=1

## method defined elsewhere

########################################################################

load("RData/data.expression.ortho.RData")

########################################################################

tpm.orthopc=exportho.mr[which(exportho.mr$GeneType=="protein_coding"), which(!colnames(exportho.mr)%in%c("ID", "GeneType"))]
tpm.ortholnc=exportho.mr[which(exportho.mr$GeneType=="lncRNA"), which(!colnames(exportho.mr)%in%c("ID", "GeneType"))]

samples=colnames(tpm.orthopc)
species=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))
tissues=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[2]))
ages=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[3]))
ages=unlist(lapply(ages, function(x) substr(x, 1, nchar(x)-1)))

tpm.orthopc=tpm.orthopc[,samples]
tpm.ortholnc=tpm.ortholnc[,samples]


########################################################################

compute.correlations <- function(expdata, tissues, ages, species, method="spearman"){
  samples.order=c()

  cor.between=c()
  cor.within=c()
  
  for(tiss in c("Brain", "Kidney", "Liver", "Testis")){
    for(age in c("EarlyEmbryo", "LateEmbryo", "Newborn", "Adult", "Aged")){
      
      samples.mouse=samples[which(tissues==tiss & ages==age & species=="Mouse")]
      samples.rat=samples[which(tissues==tiss & ages==age & species=="Rat")]
      
      if(length(samples.mouse)>=2 & length(samples.rat)>=2){

        ## select genes >=minexp in at least one sample

        samples.order=c(samples.order, paste(tiss, age, sep="_"))
        ## compute mean by species
        
        mean.mouse=apply(expdata[, samples.mouse], 1, mean)
        mean.rat=apply(expdata[, samples.rat], 1, mean)
                
        ## correlations between species

        cor.between=c(cor.between, cor(log2(mean.mouse+1), log2(mean.rat+1), use="complete.obs", method=method))
    
        ## correlations within species

        meanwithin.mouse=mean(unlist(lapply(1:(length(samples.mouse)-1), function(i) unlist(lapply((i+1):length(samples.mouse), function(j) cor(log2(expdata[,samples.mouse[i]]+1),log2(expdata[,samples.mouse[j]]+1), use="complete.obs", method=method))))))
        meanwithin.rat=mean(unlist(lapply(1:(length(samples.rat)-1), function(i) unlist(lapply((i+1):length(samples.rat), function(j) cor(log2(expdata[,samples.rat[i]]+1),log2(expdata[,samples.rat[j]]+1), use="complete.obs", method=method))))))
        
        cor.within=c(cor.within, mean(c(meanwithin.mouse, meanwithin.rat)))
        
      }
    }
  }
  
  names(cor.between)=samples.order
  names(cor.within)=samples.order

  return(list("between"=cor.between, "within"=cor.within))
}

########################################################################

corpctpm=compute.correlations(tpm.orthopc, tissues, ages, species, method=method)
corpc.between.tpm=corpctpm[["between"]]
corpc.within.tpm=corpctpm[["within"]]

corlnctpm=compute.correlations(tpm.ortholnc, tissues, ages, species, method=method)
corlnc.between.tpm=corlnctpm[["between"]]
corlnc.within.tpm=corlnctpm[["within"]]

########################################################################

nbrand=100
nbsamples=length(corpc.between.tpm)
names=names(corpc.between.tpm)

bootstrap.pc.between.tpm=matrix(rep(NA, nbrand*nbsamples), nrow=nbrand)
bootstrap.pc.within.tpm=matrix(rep(NA, nbrand*nbsamples), nrow=nbrand)

bootstrap.lnc.between.tpm=matrix(rep(NA, nbrand*nbsamples), nrow=nbrand)
bootstrap.lnc.within.tpm=matrix(rep(NA, nbrand*nbsamples), nrow=nbrand)

for(i in 1:nbrand){
  print(i)
  
  rand.pc.tpm=sample(1:dim(tpm.orthopc)[1], size=dim(tpm.orthopc)[1], replace=T)
  rand.lnc.tpm=sample(1:dim(tpm.ortholnc)[1], size=dim(tpm.ortholnc)[1], replace=T)

  rand.corpctpm=compute.correlations(tpm.orthopc[rand.pc.tpm,], tissues, ages, species, method=method)
  rand.corlnctpm=compute.correlations(tpm.ortholnc[rand.lnc.tpm,], tissues, ages, species, method=method)

  bootstrap.pc.between.tpm[i,]=rand.corpctpm[["between"]]
  bootstrap.pc.within.tpm[i,]=rand.corpctpm[["within"]]

  bootstrap.lnc.between.tpm[i,]=rand.corlnctpm[["between"]]
  bootstrap.lnc.within.tpm[i,]=rand.corlnctpm[["within"]]
}

colnames(bootstrap.pc.between.tpm)=names
colnames(bootstrap.pc.within.tpm)=names

colnames(bootstrap.lnc.between.tpm)=names
colnames(bootstrap.lnc.within.tpm)=names

########################################################################

save(list=c("corpc.between.tpm", "corlnc.between.tpm", "corpc.within.tpm", "corlnc.within.tpm", "bootstrap.pc.between.tpm", "bootstrap.pc.within.tpm", "bootstrap.lnc.between.tpm", "bootstrap.lnc.within.tpm"), file=paste("RData/data.expression.conservation.", method,".RData",sep=""))

########################################################################
########################################################################

