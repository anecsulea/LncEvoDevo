########################################################################

path="LncEvoDevo/"
pathExpression=paste(path, "results/expression_estimation/", sep="")

########################################################################

for(sp in c("Mouse", "Rat", "Chicken")){

  print(sp)

  exp.allregions=read.table(paste(pathExpression, sp, "/AllSamples_UniqueReadCounts_StringTie_MainStrain.txt", sep=""), h=T, stringsAsFactors=F)
  exp.noov=read.table(paste(pathExpression, sp, "/AllSamples_UniqueReadCounts_NonOverlappingExonBlocks_StringTie_MainStrain.txt", sep=""), h=T, stringsAsFactors=F)

  rownames(exp.allregions)=exp.allregions$GeneID
  rownames(exp.noov)=exp.noov$GeneID
  
  info.allregions=exp.allregions[,c("GeneID","ExonicLength")]
  info.noov=exp.noov[,c("GeneID","ExonicLength")]

  exp.allregions=exp.allregions[,setdiff(colnames(exp.allregions), c("GeneID","ExonicLength"))]
  exp.noov=exp.noov[,setdiff(colnames(exp.noov), c("GeneID","ExonicLength"))]
  
  samples=colnames(exp.allregions)
  exp.noov=exp.noov[,samples]
  
  exp.allregions=as.matrix(exp.allregions)
  exp.noov=as.matrix(exp.noov)

  total.allregions=apply(exp.allregions,1,sum)
  total.noov=apply(exp.noov, 1, sum)

  names(total.allregions)=rownames(exp.allregions)
  names(total.noov)=rownames(exp.noov)

  results=data.frame("GeneID"=info.allregions$GeneID, "ExonicLength"=info.allregions$ExonicLength, stringsAsFactors=F)
  rownames(results)=results$GeneID
  
  results$ExonicLengthNoOverlap=rep(0, dim(results)[1])
  results[info.noov$GeneID,"ExonicLengthNoOverlap"]=info.noov$ExonicLength

  results$ExpressionCorrelation=rep(NA, dim(results)[1])

  results[info.noov$GeneID, "ExpressionCorrelation"]=unlist(lapply(info.noov$GeneID, function(x) cor(as.numeric(exp.allregions[x,]), as.numeric(exp.noov[x,]), method="spearman")))

  results$TotalReadCount=total.allregions[results$GeneID]
  results$ReadCountNoOverlap=rep(0, dim(results)[1])

  results[info.noov$GeneID, "ReadCountNoOverlap"]=total.noov[info.noov$GeneID]
          
  write.table(results, file=paste(pathExpression, sp, "/Comparison_AllRegions_NonOverlappingExonBlocks_StringTie_MainStrain.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  
}

########################################################################

