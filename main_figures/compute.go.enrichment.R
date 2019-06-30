##################################################################

compute.go.enrichment <- function(targetlist, background, categories, genelist){
  pvalues=c()
  enrichment=c()

  background=unique(c(background, targetlist))

  N=length(background)
  n=length(targetlist)

  allb=c()
  allB=c()
  
  for(cat in categories[,1]){
    B=length(intersect(background, genelist[[cat]]))
    b=length(intersect(targetlist, genelist[[cat]]))

    allb=c(allb, b)
    allB=c(allB, B)

    if(b>1){
      pval=1-phyper(q=b-1, k=n, m=B, n=N-B)
    } else{
      pval=1
    }

    e=((b/n)/(B/N))

    pvalues=c(pvalues, pval)
    enrichment=c(enrichment, e)
  }

  results=data.frame("GOAccession"=categories[,1], "GOName"=categories[,2],  "PValue"=pvalues, "Enrichment"=enrichment, "N"=rep(N, dim(categories)[1]), "n"=rep(n, dim(categories)[1]), "b"=allb, "B"=allB, stringsAsFactors=F)
  results$FDR=p.adjust(results$PValue,method="BH")

  return(results)

}

##################################################################
