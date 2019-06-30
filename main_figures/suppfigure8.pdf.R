###################################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE
  
  minTPM=1
  maxFDR=0.01
}


###################################################################################

if(load==TRUE){
  for(sp in c("Mouse", "Rat", "Chicken")){
    load(paste("RData/data.annotations.", sp, ".RData", sep=""))
    assign(paste("pc", tolower(sp), sep="."), pc)
    assign(paste("lnc", tolower(sp), sep="."), lnc)
    assign(paste("allinfo", tolower(sp), sep="."), allinfo)
    
    load(paste("RData/data.expression.",sp,".RData", sep=""))
    assign(paste("normtpm", tolower(sp), sep="."), normtpm)

    load(paste("RData/data.tpm.stats.",sp,".RData", sep=""))
    assign(paste("stats", tolower(sp), sep="."), stats) 

    
  }

  load=FALSE
}

###################################################################################

if(prepare==TRUE){

  alldata=list()
  
  for(sp in c("Mouse", "Rat", "Chicken")){

    alldata[[sp]]=list()
    
    this.normtpm=get(paste("normtpm", tolower(sp), sep="."))
    
    this.info=get(paste("allinfo", tolower(sp), sep="."))
    
    this.pc=get(paste("pc", tolower(sp), sep="."))
    this.lnc=get(paste("lnc", tolower(sp), sep="."))
    
    this.spliced=this.info$GeneID[which(this.info$NbExons>1)]
    this.bidirpc=this.info$GeneID[which(!is.na(this.info$BidirectionalPromoterProteinCoding))]
    this.bidirother=setdiff(this.info$GeneID[which(!is.na(this.info$BidirectionalPromoter))], this.bidirpc)
    this.enhancers=this.info$GeneID[which(this.info$OverlapEncodeEnhancer=="Yes")]
    this.antisense=this.info$GeneID[which(this.info$RegionPCAntisense!="intergenic")]
    
    maxexp=log2(apply(this.normtpm, 1, function(x) max(x, na.rm=T))+1)
    names(maxexp)=rownames(this.normtpm)
    
    alldata[[sp]][["pc"]]=this.pc
    alldata[[sp]][["lnc"]]=this.lnc
    alldata[[sp]][["spliced"]]=this.spliced
    alldata[[sp]][["bidirpc"]]=this.bidirpc
    alldata[[sp]][["bidirother"]]=this.bidirother
    alldata[[sp]][["enhancers"]]=this.enhancers
    alldata[[sp]][["antisense"]]=this.antisense
    alldata[[sp]][["maxexp"]]=maxexp
    
  }

  prepare=FALSE
}

###################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure8.pdf", sep=""), width=6.0, height=8.5)

###################################################################################

## layout

m=matrix(rep(NA,3*9), nrow=3)

m[1,]=c(rep(1, 2), rep(2, 7))

m[2,]=c(rep(3, 2), rep(4, 7))

m[3,]=c(rep(5, 2), rep(6, 7))

layout(m)

###################################################################################

## maximum expression level

labels=toupper(letters[1:6])
indexplot=0

for(sp in c("Mouse", "Rat", "Chicken")){

  pc=alldata[[sp]][["pc"]]
  lnc=alldata[[sp]][["lnc"]]
  spliced=alldata[[sp]][["spliced"]]
  bidirpc=alldata[[sp]][["bidirpc"]]
  bidirother=alldata[[sp]][["bidirother"]]
  enhancers=alldata[[sp]][["enhancers"]]
  antisense=alldata[[sp]][["antisense"]]
  maxexp=alldata[[sp]][["maxexp"]]
  
  ## max expression level, pc against lnc
  indexplot=indexplot+1

  ylim=c(0, 13)
  xlim=c(0.5, 2.5)

  par(mar=c(7.75, 2.75, 1.5, 0.1))
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)
  
  boxplot(maxexp[pc],at=1, col="indianred", axes=F, add=T, outline=F, boxwex=0.7, notch=T)
  boxplot(maxexp[lnc],at=2, col="steelblue", axes=F, add=T, outline=F, boxwex=0.7, notch=T)

  axis(side=1, at=c(1:2), mgp=c(3, 0.5, 0), labels=rep("",2))
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.9)
  
  mtext("maximum expression (log2 TPM)", side=2, cex=0.65, line=1.75)

  mtext(c("protein\ncoding", "lncRNA"), side=1, at=1:2, line=0.75, cex=0.6, las=2)

  mtext(labels[indexplot], side=3, at=-0.22, line=0.5, font=2, cex=0.95)
  
  ## max expression level, lnc types

  indexplot=indexplot+1
  
  ylim=c(0, 8)
  xlim=c(0.5, 11.5)

  xpos=c(1, 2.25, 3.25, 4.5,5.5,6.5,  7.75,8.75, 10,11)
  
  par(mar=c(7.75, 2.75, 1.5, 0.1))
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

  boxplot(maxexp[lnc], at=xpos[1], col="steelblue", axes=F, add=T, outline=F, boxwex=0.7, notch=T)

  boxplot(maxexp[intersect(lnc, spliced)], at=xpos[2], col="steelblue", axes=F, add=T, outline=F, boxwex=0.7, notch=T)
  boxplot(maxexp[setdiff(lnc, spliced)], at=xpos[3], col="steelblue", axes=F, add=T, outline=F, boxwex=0.7, notch=T)

  boxplot(maxexp[intersect(lnc, bidirpc)], at=xpos[4], col="steelblue", axes=F, add=T, outline=F, boxwex=0.7, notch=T)
  boxplot(maxexp[intersect(lnc, bidirother)], at=xpos[5], col="steelblue", axes=F, add=T, outline=F, boxwex=0.7, notch=T)
  boxplot(maxexp[setdiff(lnc, c(bidirpc, bidirother))], at=xpos[6], col="steelblue", axes=F, add=T, outline=F, boxwex=0.7, notch=T)
  
  boxplot(maxexp[intersect(lnc, antisense)], at=xpos[7], col="steelblue", axes=F, add=T, outline=F, boxwex=0.7, notch=T)
  boxplot(maxexp[setdiff(lnc, antisense)], at=xpos[8], col="steelblue", axes=F, add=T, outline=F, boxwex=0.7, notch=T)

  if(length(enhancers)>0){
    boxplot(maxexp[intersect(lnc, enhancers)], at=xpos[9], col="steelblue", axes=F, add=T, outline=F, boxwex=0.7, notch=T)
    boxplot(maxexp[setdiff(lnc, enhancers)], at=xpos[10], col="steelblue", axes=F, add=T, outline=F, boxwex=0.7, notch=T)
  }


  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.9)

  if(length(enhancers)>0){
    axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("",length(xpos)))
    mtext(c("all", "spliced", "unspliced", "bidirectional promoter\nprotein-coding","bidirectional\npromoter, other", "unidirectional\npromoter", "antisense", "intergenic", "enhancer\noverlap", "no enhancer\noverlap"), at=xpos, side=1, cex=0.6, line=0.5, las=2)
    mtext(tolower(sp), side=3, cex=0.65, line=-0.2, at=max(xpos))
  } else{
    axis(side=1, at=xpos[1:(length(xpos)-2)], mgp=c(3, 0.5, 0), labels=rep("",length(xpos)-2))
    mtext(c("all", "spliced", "unspliced", "bidirectional promoter\n protein-coding","bidirectional\npromoter, other", "unidirectional\npromoter", "antisense", "intergenic"), at=xpos[1:8], side=1, cex=0.6, line=0.5, las=2)
    mtext(tolower(sp), side=3, cex=0.65, line=-0.2, at=max(xpos[1:(length(xpos)-2)]))
  }

  mtext("maximum expression (log2 TPM)", side=2, cex=0.65, line=1.75)

  
  
  
  mtext(labels[indexplot], side=3, at=-0.7, line=0.5, font=2, cex=0.95)
}



###################################################################################


dev.off()

###################################################################################
