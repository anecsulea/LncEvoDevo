

normalization <- function(expdata,nbgenes=100, method="HK"){
  if(method=="HK"){
    
  ### get genes that are expressed in all samples
    
  expdata.nonzero=expdata[which(apply(expdata,1,min)>0 ),]
  
  ## transform the expression levels into ranks
  
  expdata.ranks=apply(expdata.nonzero,2,rank)
  
  ## compute the variance of the ranks for each gene
  
  expdata.ranks.var=apply(expdata.ranks,1,var,na.rm=T)
  
  ## rank the genes according to their variance
  
  expdata.nonzero.consrank=rank(expdata.ranks.var)
  
  ## compute the median rank over all samples, for each gene
  
  median.rank=apply(expdata.ranks,1,median,na.rm=T)
  
  ## get the genes which have a median rank in the 25%-75% range
  
  interquartile=median.rank>(0.25*length(median.rank)) & median.rank<(0.75*length(median.rank))
  
  ## get the house-keeping genes: the n genes with the most conserved ranks, among those that fall in the interquartile range
  
  hkgenes=names(sort(expdata.nonzero.consrank[interquartile])[1:nbgenes])
  
  ## compute the normalization coefficient for each sample =  the median of the RPKM of the hkgenes in that sample
  
  normcoeff=apply(expdata.nonzero[hkgenes,],2,median,na.rm=T)
  
  ## we want to bring all of the medians at the average of the medians (computed on all expressed genes)
  
 ##  normcoeff=normcoeff/mean(apply(expdata.nonzero,2,median))

  normcoeff=normcoeff/mean(normcoeff)
  
  ## finally, normalize the data
  
  expdata.norm=t(t(expdata)/normcoeff)

  expdata.norm=as.data.frame(expdata.norm)

  results=list("expdata.norm"=expdata.norm, "hk.genes"=hkgenes, "normcoeff"=normcoeff)
  return(results)
  
}
  else{
    print("other methods are not implemented yet")
  }
  
}
