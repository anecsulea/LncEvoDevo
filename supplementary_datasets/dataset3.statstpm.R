##########################################################################

path="LncEvoDevo/"
pathExpression=paste(path, "results/expression_estimation/", sep="")
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathResults=paste(path, "supplementary_datasets/SupplementaryDataset3/", sep="")

options(scipen=999) ## remove scientific notation ## options(scipen=0) to get it back

release=94

### average TPM values per tissue and age, tau values

###############################################################################

compute.tau <- function(exp){
  if(max(exp)==0){
    return(NA)
  }
  
  n=length(exp)
  newexp=exp/max(exp)

  tau=sum(1-newexp)/(n-1)

  return(tau)
}

##########################################################################

for(sp in c("Mouse", "Rat", "Chicken")){
  print(sp)

  
  norm.tpm=read.table(paste(pathExpression, sp, "/AllSamples_KallistoNormalizedTPM_StringTie_MainStrain.txt", sep=""), h=T, stringsAsFactors=F)
  rownames(norm.tpm)=norm.tpm$GeneID

  norm.tpm=norm.tpm[,-which(colnames(norm.tpm)%in%c("GeneID"))]

  samples=colnames(norm.tpm)
  tissue=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))
  age=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[2]))
  age=unlist(lapply(age, function(x) substr(x, 1, nchar(x)-1)))
  tissage=as.factor(paste(tissue, age, sep="_"))

  
  tpm.tissage=t(apply(norm.tpm, 1, function(x) tapply(as.numeric(x), tissage, mean)))
  
  maxexp=apply(tpm.tissage, 1, max)
  maxsample=apply(tpm.tissage, 1, function(x) colnames(tpm.tissage)[which.max(x)])
  maxsample[which(maxexp==0)]=NA
  
  colnames(tpm.tissage)=paste("MeanTPM", colnames(tpm.tissage), sep=".")

  tau.tpm=apply(tpm.tissage, 1, function(x) compute.tau(x))

  ## compute tissue/stage specificity with just 4 developmental stages, combining young adult and aged adult

  age4=age
  age4[which(age4%in%c("Adult", "Aged"))]="AdultAged"
  tissage4=as.factor(paste(tissue, age4, sep="_"))

  tpm.tissage4=t(apply(norm.tpm, 1, function(x) tapply(as.numeric(x), tissage4, mean)))
  maxexp4=apply(tpm.tissage4, 1, max)
  maxsample4=apply(tpm.tissage4, 1, function(x) colnames(tpm.tissage4)[which.max(x)])
  maxsample4[which(maxexp4==0)]=NA
  
  colnames(tpm.tissage4)=paste("MeanTPM", colnames(tpm.tissage4), sep=".")

  tau.tpm4=apply(tpm.tissage4, 1, function(x) compute.tau(x))
  
  results=data.frame("GeneID"=rownames(norm.tpm), tpm.tissage, tpm.tissage4[,setdiff(colnames(tpm.tissage4), colnames(tpm.tissage))], "ExpressionSpecificity"=tau.tpm, "MaxSample"=maxsample,  "MaxExpression"=maxexp, "ExpressionSpecificity4Stages"=tau.tpm4, "MaxSample4Stages"=maxsample4, "MaxExpression4Stages"=maxexp4, stringsAsFactors=F)
  
  write.table(results, file=paste(pathResults, "Statistics_AverageTPM_",sp,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

}

##########################################################################
