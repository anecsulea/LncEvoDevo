###################################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE
  
  minTPM=1
}


###################################################################################

if(load==TRUE){
  sp="Mouse"
  phasttype="placental"
  
  load(paste("RData/data.annotations.", sp, ".RData", sep=""))
 
  load(paste("RData/data.diffexp.",sp,".RData", sep=""))
  de.global=de.global.allreads

  load(paste("RData/data.phastcons.", sp,".", phasttype,".RData", sep=""))

  load(paste("RData/data.tpm.stats.", sp,".RData", sep=""))
  avgtpm=stats

  load=FALSE
}

###################################################################################

if(prepare==TRUE){
  
  spliced=allinfo$GeneID[which(allinfo$NbExons>1)]
  bidirpc=allinfo$GeneID[which(!is.na(allinfo$BidirectionalPromoterProteinCoding))]
  bidirother=setdiff(allinfo$GeneID[which(!is.na(allinfo$BidirectionalPromoter))], bidirpc)
  unidir=setdiff(allinfo$GeneID, c(bidirpc, bidirother))
  enhancers=allinfo$GeneID[which(allinfo$OverlapEncodeEnhancer=="Yes")]
  antisense=allinfo$GeneID[which(allinfo$RegionPCAntisense!="intergenic")]
  intergenic=allinfo$GeneID[which(allinfo$RegionPCAntisense=="intergenic")]


  exons.score=phastcons[["exons"]]
  introns.score=phastcons[["introns"]]
  promoters.score=phastcons[["promoters"]]
  splice5.score=phastcons[["splice5"]]
  splice3.score=phastcons[["splice3"]]
  splicesites.score=phastcons[["splicesites"]]
  intergenic.score=phastcons[["intergenic"]]
  mean.intergenic=phastcons[["mean.intergenic"]]
    
  
  ## take only genes above noise expression level, in at least one sample

  
  maxsample=avgtpm[,"MaxSample"]
  maxtissue=unlist(lapply(maxsample, function(x) unlist(strsplit(x, split="_"))[1]))
  maxage=unlist(lapply(maxsample, function(x) unlist(strsplit(x, split="_"))[2]))

  names(maxsample)=rownames(avgtpm)
  names(maxtissue)=rownames(avgtpm)
  names(maxage)=rownames(avgtpm)
  
  ## phastcons scores
 
  medians.promoters.unidir=c()
  medians.promoters.bidirpc=c()
  medians.promoters.bidirother=c()
  
  ci.low.promoters.unidir=c()
  ci.low.promoters.bidirpc=c()
  ci.low.promoters.bidirother=c()

  ci.high.promoters.unidir=c()
  ci.high.promoters.bidirpc=c()
  ci.high.promoters.bidirother=c()

  
  ## median values per age
  
  diff.exons.splicesites=list()
  diff.exons.promoters=list()

  stages=c()
  
  for(tiss in tissue.order){
    for(age in age.order){
      this.pc=pc[which(avgtpm[pc,paste("MeanTPM.", tiss, "_",age,sep="")]>=minTPM)]
      this.lnc=lnc[which(avgtpm[lnc,paste("MeanTPM.", tiss, "_",age,sep="")]>=minTPM)]

      if(length(this.pc)>0 & length(this.lnc)>0){
        diff.exons.splicesites[[paste(tiss, age, sep="_")]]=list("pc"=(exons.score[this.pc]-splicesites.score[this.pc]), "lnc"=(exons.score[this.lnc]-splicesites.score[this.lnc]))
        
        diff.exons.promoters[[paste(tiss, age, sep="_")]]=list("pc"=(exons.score[this.pc]-promoters.score[this.pc]), "lnc"=(exons.score[this.lnc]-promoters.score[this.lnc]))
        
        this.lnc.unidir=intersect(this.lnc, unidir)
        this.lnc.bidirother=intersect(this.lnc, bidirother)
        this.lnc.bidirpc=intersect(this.lnc, bidirpc)
        
        conf.lnc.unidir=boxplot(promoters.score[this.lnc.unidir], plot=F)$conf
        conf.lnc.bidirother=boxplot(promoters.score[this.lnc.bidirother], plot=F)$conf
        conf.lnc.bidirpc=boxplot(promoters.score[this.lnc.bidirpc], plot=F)$conf
        
        medians.promoters.unidir=c(medians.promoters.unidir, median(promoters.score[this.lnc.unidir], na.rm=T))
        medians.promoters.bidirpc=c(medians.promoters.bidirpc, median(promoters.score[this.lnc.bidirpc], na.rm=T))
        medians.promoters.bidirother=c(medians.promoters.bidirother, median(promoters.score[this.lnc.bidirother], na.rm=T))

        ci.low.promoters.unidir=c(ci.low.promoters.unidir, conf.lnc.unidir[1])
        ci.high.promoters.unidir=c(ci.high.promoters.unidir, conf.lnc.unidir[2])

        ci.low.promoters.bidirpc=c(ci.low.promoters.bidirpc, conf.lnc.bidirpc[1])
        ci.high.promoters.bidirpc=c(ci.high.promoters.bidirpc, conf.lnc.bidirpc[2])

        ci.low.promoters.bidirother=c(ci.low.promoters.bidirother, conf.lnc.bidirother[1])
        ci.high.promoters.bidirother=c(ci.high.promoters.bidirother, conf.lnc.bidirother[2])
        
        stages=c(stages, paste(tiss, age, sep="_"))
        
      }
    }
  }

  names(medians.promoters.unidir)=stages
  names(medians.promoters.bidirpc)=stages
  names(medians.promoters.bidirother)=stages
  
  names(ci.low.promoters.unidir)=stages
  names(ci.low.promoters.bidirpc)=stages
  names(ci.low.promoters.bidirother)=stages
  
  names(ci.high.promoters.unidir)=stages
  names(ci.high.promoters.bidirpc)=stages
  names(ci.high.promoters.bidirother)=stages
  
  prepare=FALSE
}

###################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure11.pdf", sep=""), width=6.85, height=5.75)

m=matrix(rep(NA,20*12), nrow=20)

for(i in 1:10){
  m[i,]=c(rep(1,4), rep(2,4), rep(3,4))
}

for(i in 11:20){ 
  m[i,]=c(rep(4,6), rep(5,6))
}

layout(m)

###################################################################################

samples=kronecker(tissue.order, age.order, paste, sep="_")
xpos=1:20
xpos=setdiff(xpos, 16)
xlim=c(0,20.5)

###################################################################################

## promoters, lnc with unidirectional promoters

ylim=range(c(ci.low.promoters.unidir, ci.high.promoters.unidir))
ylim[1]=0.05
ylim[2]=0.27

par(mar=c(3.5,2.75,2.1,0.5))
plot(1, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F)

abline(h=mean.intergenic,  col="black",lty=3)

this.tissues=unlist(lapply(names(medians.promoters.unidir), function(x) unlist(strsplit(x, split="_"))[1]))

segments(xpos, ci.low.promoters.unidir, xpos, ci.high.promoters.unidir, col=col.tissues[this.tissues])
points(xpos, medians.promoters.unidir, pch=20, col=col.tissues[this.tissues])

axis(side=2,  mgp=c(3,0.5,0), cex.axis=0.85)
axis(side=1, at=1:20, labels=rep("",20), mgp=c(3,0.75,0), cex.axis=0.85)
mtext(rep(1:5, 4), side=1, at=1:20, cex=0.6, line=0.5)
mtext("developmental stage", side=1, line=1.5, cex=0.65)
mtext("sequence conservation", side=2, line=2, cex=0.65)

mtext("unidirectional promoters", side=3, line=0, cex=0.65)
mtext("A", side=3, at=-4.35, line=0.5, font=2, cex=0.9)

###################################################################################

## promoters, lnc with bidirectional promoters, protein-coding gene 

ylim=range(c(ci.low.promoters.bidirpc, ci.high.promoters.bidirpc))
ylim[1]=0.05
ylim[2]=0.35

par(mar=c(3.5,2.75,2.1,0.5))
plot(1, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F)

abline(h=mean.intergenic,  col="black",lty=3)

this.tissues=unlist(lapply(names(medians.promoters.bidirpc), function(x) unlist(strsplit(x, split="_"))[1]))

segments(xpos, ci.low.promoters.bidirpc, xpos, ci.high.promoters.bidirpc, col=col.tissues[this.tissues])
points(xpos, medians.promoters.bidirpc, pch=20, col=col.tissues[this.tissues])

axis(side=2,  mgp=c(3,0.5,0), cex.axis=0.85)
axis(side=1, at=1:20, labels=rep("",20), mgp=c(3,0.75,0), cex.axis=0.85)
mtext(rep(1:5, 4), side=1, at=1:20, cex=0.6, line=0.5)
mtext("developmental stage", side=1, line=1.5, cex=0.65)
mtext("sequence conservation", side=2, line=2, cex=0.65)

mtext("bidirectional promoters, protein-coding", side=3, line=0, cex=0.65)
mtext("B", side=3, at=-4.35, line=0.5, font=2, cex=0.9)

###################################################################################

## promoters, lnc with bidirectional promoters, other gene

ylim=range(c(ci.low.promoters.bidirother, ci.high.promoters.bidirother))
ylim[1]=0.05
ylim[2]=0.27


par(mar=c(3.5,2.75,2.1,0.5))
plot(1, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F)

abline(h=mean.intergenic,  col="black",lty=3)

this.tissues=unlist(lapply(names(medians.promoters.bidirother), function(x) unlist(strsplit(x, split="_"))[1]))

segments(xpos, ci.low.promoters.bidirother, xpos, ci.high.promoters.bidirother, col=col.tissues[this.tissues])
points(xpos, medians.promoters.bidirother, pch=20, col=col.tissues[this.tissues])

axis(side=2,  mgp=c(3,0.5,0), cex.axis=0.85)
axis(side=1, at=1:20, labels=rep("",20), mgp=c(3,0.75,0), cex.axis=0.85)
mtext(rep(1:5, 4), side=1, at=1:20, cex=0.6, line=0.5)
mtext("developmental stage", side=1, line=1.5, cex=0.65)
mtext("sequence conservation", side=2, line=2, cex=0.65)


mtext("bidirectional promoters, other", side=3, line=0, cex=0.65)
mtext("C", side=3, at=-4.35, line=0.5, font=2, cex=0.9)

###################################################################################

## difference between exons and promoters

samples=kronecker(tissue.order, age.order, paste, sep="_")
xpos=1:20
names(xpos)=samples

ylim=c(-0.11,0.11)

par(mar=c(3.5,2.5,2.1,1.5))

plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(0.5, 20.5), ylim=ylim)

abline(h=0, lty=2, col="gray40")

for(tiss in tissue.order){
  for(age in age.order){
    this.sample=paste(tiss, age, sep="_")

    if(this.sample%in%names(diff.exons.promoters)){
      this.diff.pc=diff.exons.promoters[[paste(tiss, age, sep="_")]][["pc"]]
      this.diff.lnc=diff.exons.promoters[[paste(tiss, age, sep="_")]][["lnc"]]
      
      this.x=xpos[paste(tiss, age, sep="_")]
      
      median.lnc=median(this.diff.lnc, na.rm=T)
      ci.lnc=boxplot(this.diff.lnc, plot=F)$conf
      
      segments(this.x, ci.lnc[1], this.x, ci.lnc[2], col=col.tissues[tiss])
      points(this.x, median.lnc, col=col.tissues[tiss], pch=20)
    }
  }
}


axis(side=2,  mgp=c(3,0.5,0), cex.axis=0.7)
axis(side=1, at=1:20, labels=rep("",20), mgp=c(3,0.75,0), cex.axis=0.7)

mtext(rep(1:5, 4), at=1:20, cex=0.55, line=0.5, side=1)
mtext("sequence conservation difference", side=2, line=1.5, cex=0.65)
mtext("lncRNA exons - promoters", side=3, line=0, cex=0.65)
mtext("developmental stage", side=1, line=1.5, cex=0.65)

legend("topleft", legend=shortname.tiss, border=col.tissues, fill=col.tissues, horiz=F, bty="n", xpd=NA, inset=c(0, 0.05), cex=0.95)

mtext("D", side=3, at=-2.25, line=0.5, font=2, cex=0.9)

###################################################################################


## difference between exons and splicesites

samples=kronecker(tissue.order, age.order, paste, sep="_")
xpos=1:20
names(xpos)=samples

ylim=c(-0.11,0.11)

par(mar=c(3.5,2.5,2.1,1.5))

plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(0.5, 20.5), ylim=ylim)

abline(h=0, lty=2, col="gray40")

for(tiss in tissue.order){
  for(age in age.order){
    this.sample=paste(tiss, age, sep="_")

    if(this.sample%in%names(diff.exons.splicesites)){
      this.diff.pc=diff.exons.splicesites[[paste(tiss, age, sep="_")]][["pc"]]
      this.diff.lnc=diff.exons.splicesites[[paste(tiss, age, sep="_")]][["lnc"]]
      
      this.x=xpos[paste(tiss, age, sep="_")]
      
      median.lnc=median(this.diff.lnc, na.rm=T)
      ci.lnc=boxplot(this.diff.lnc, plot=F)$conf
      
      segments(this.x, ci.lnc[1], this.x, ci.lnc[2], col=col.tissues[tiss])
      points(this.x, median.lnc, col=col.tissues[tiss], pch=20)
    }
  }
}


axis(side=2,  mgp=c(3,0.5,0), cex.axis=0.7)
axis(side=1, at=1:20, labels=rep("",20), mgp=c(3,0.75,0), cex.axis=0.7)

mtext(rep(1:5, 4), at=1:20, cex=0.55, line=0.5, side=1)
mtext("sequence conservation difference", side=2, line=1.5, cex=0.65)
mtext("lncRNA exons - splice sites", side=3, line=0, cex=0.65)
mtext("developmental stage", side=1, line=1.5, cex=0.65)


mtext("E", side=3, at=-2.25, line=0.5, font=2, cex=0.9)

###################################################################################

dev.off()

###################################################################################
