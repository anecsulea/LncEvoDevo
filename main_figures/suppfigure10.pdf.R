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
  
  allinfo=allinfo[which(allinfo$Chr%in%c(as.character(1:19), "X", "Y")),]
  pc=intersect(pc, allinfo$GeneID)
  lnc=intersect(lnc, allinfo$GeneID)
  
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

  ## expression patterns
  
  maxsample=avgtpm[,"MaxSample"]
  maxtissue=unlist(lapply(maxsample, function(x) unlist(strsplit(x, split="_"))[1]))
  maxage=unlist(lapply(maxsample, function(x) unlist(strsplit(x, split="_"))[2]))

  names(maxsample)=rownames(avgtpm)
  names(maxtissue)=rownames(avgtpm)
  names(maxage)=rownames(avgtpm)


  ## median values per age
  
  allvalues=list()
  
  for(type in c("exons", "splicesites", "promoters")){
        
    this.cons=get(paste(type, ".score", sep=""))
  
    allvalues[[type]]=list()
    
    medians.pc=c()
    medians.lnc=c()
    ci.low.pc=c()
    ci.high.pc=c()
    ci.low.lnc=c()
    ci.high.lnc=c()
    
    nb.pc=c()
    nb.lnc=c()

    stages=c()
  
    
    for(tiss in tissue.order){
      for(age in age.order){
        this.pc=pc[which(avgtpm[pc,paste("MeanTPM.", tiss, "_",age,sep="")]>=minTPM)]
        this.lnc=lnc[which(avgtpm[lnc,paste("MeanTPM.", tiss, "_",age,sep="")]>=minTPM)]
     
        this.pc=intersect(this.pc, names(this.cons))
        this.lnc=intersect(this.lnc, names(this.cons))

        if(length(this.pc)>0 & length(this.lnc)>0){
          medians.pc=c(medians.pc, median(this.cons[this.pc], na.rm=T))
          medians.lnc=c(medians.lnc, median(this.cons[this.lnc], na.rm=T))

          c.pc=boxplot(this.cons[this.pc],plot=F)$conf
          ci.low.pc=c(ci.low.pc, c.pc[1,1])
          ci.high.pc=c(ci.high.pc, c.pc[2,1])

          c.lnc=boxplot(this.cons[this.lnc],plot=F)$conf
          ci.low.lnc=c(ci.low.lnc, c.lnc[1,1])
          ci.high.lnc=c(ci.high.lnc, c.lnc[2,1])
          
          
          nb.pc=c(nb.pc, length(this.pc))
          nb.lnc=c(nb.lnc, length(this.lnc))

          stages=c(stages, paste(tiss, age, sep="_"))
       
        }
      }

      names(medians.pc)=stages
      names(medians.lnc)=stages
      
      names(ci.low.pc)=stages
      names(ci.low.lnc)=stages
      
      names(ci.high.pc)=stages
      names(ci.high.lnc)=stages

      names(nb.pc)=stages
      names(nb.lnc)=stages

      allvalues[[type]][["medians.pc"]]=medians.pc
      allvalues[[type]][["medians.lnc"]]=medians.lnc

      allvalues[[type]][["ci.low.pc"]]=ci.low.pc
      allvalues[[type]][["ci.low.lnc"]]=ci.low.lnc
      
      allvalues[[type]][["ci.high.pc"]]=ci.high.pc
      allvalues[[type]][["ci.high.lnc"]]=ci.high.lnc

      allvalues[[type]][["nb.pc"]]=nb.pc
      allvalues[[type]][["nb.lnc"]]=nb.lnc

      median.ig=median(intergenic.score)
      ci.ig=boxplot(intergenic.score)$conf
      ci.low.ig=ci.ig[1]
      ci.high.ig=ci.ig[2]
      

    }
  }  

  prepare=FALSE
}

###################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in


###################################################################################

## plot by types

pdf(file=paste(pathFigures, "SupplementaryFigure10.pdf", sep=""), width=6.85, height=3)

###################################################################################

## layout

m=matrix(rep(NA,1*30), nrow=1)

m[1,]=c(rep(1,6), rep(2,18), rep(3,6))

layout(m)

###################################################################################

## phastcons for exons

ylim=c(0,1)
xlim=c(0.5,3.5)

par(mar=c(5.1,3.1,2.2,0.25))

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

boxplot(exons.score[pc], boxwex=0.75, col="indianred", at=1, add=T, axes=F, outline=F, notch=T)
boxplot(exons.score[lnc], boxwex=0.75, col="steelblue", at=2, add=T, axes=F, outline=F, notch=T)
boxplot(intergenic.score, boxwex=0.75, col="gray60", at=3, add=T, axes=F, outline=F, notch=T)

axis(side=1, at=c(1,2,3), labels=rep("",3))
axis(side=2, cex.axis=0.85, mgp=c(3,0.5,0))
mtext("sequence conservation", side=2, line=2, cex=0.65)
mtext(c("protein\ncoding", "lncRNA", "intergenic"), side=1, at=1:3, line=0.75, cex=0.6, las=2)


n=c(length(pc), length(lnc))

mtext("N =", side=1, at=0, line=3.9, cex=0.6)
mtext(n, side=1, cex=0.6, line=3.9, at=c(1,2))

mtext("exons", side=3, line=0.15, cex=0.65)

mtext("A", side=3, at=-0.7, font=2, cex=0.9, line=0.75)

###################################################################################

## promoter regions


ylim=c(0,1)
xlim=c(0.6,11.5)

xpos=c(1,2, 3.5,4.5, 6,7, 8.5,9.5, 11)

par(mar=c(5.1,3.1,2.2,0.5))

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

boxplot(promoters.score[intersect(pc, unidir)], boxwex=0.7, col="indianred", at=xpos[1], add=T, axes=F, outline=F, notch=T)
boxplot(promoters.score[intersect(lnc, unidir)], boxwex=0.7, col="steelblue", at=xpos[2], add=T, axes=F, outline=F, notch=T)

boxplot(promoters.score[intersect(pc, bidirpc)], boxwex=0.7, col="indianred", at=xpos[3], add=T, axes=F, outline=F, notch=T)
boxplot(promoters.score[intersect(lnc, bidirpc)], boxwex=0.7, col="steelblue", at=xpos[4], add=T, axes=F, outline=F, notch=T)

boxplot(promoters.score[intersect(pc, bidirother)], boxwex=0.7, col="indianred", at=xpos[5], add=T, axes=F, outline=F, notch=T)
boxplot(promoters.score[intersect(lnc, bidirother)], boxwex=0.7, col="steelblue", at=xpos[6], add=T, axes=F, outline=F, notch=T)

boxplot(promoters.score[intersect(pc, enhancers)], boxwex=0.7, col="indianred", at=xpos[7], add=T, axes=F, outline=F, notch=T)
boxplot(promoters.score[intersect(lnc, enhancers)], boxwex=0.7, col="steelblue", at=xpos[8], add=T, axes=F, outline=F, notch=T)


boxplot(intergenic.score, boxwex=0.7, col="gray60", at=xpos[9], add=T, axes=F, outline=F, notch=T)


n=c(length(intersect(pc, unidir)), length(intersect(lnc, unidir)), length(intersect(pc, bidirpc)), length(intersect(lnc, bidirpc)),  length(intersect(pc, bidirother)), length(intersect(lnc, bidirother)),  length(intersect(pc, enhancers)), length(intersect(lnc, enhancers)))

mtext(n, side=1, cex=0.6, line=3.9, at=xpos)

axis(side=1, at=xpos, labels=rep("",length(xpos)))
axis(side=2, cex.axis=0.85, mgp=c(3,0.5,0))
mtext("sequence conservation", side=2, line=2, cex=0.6)

mtext("unidirectional", side=1, at=mean(xpos[1:2]), line=1.25, cex=0.6)
mtext("bidirectional,\nprotein-coding", side=1, at=mean(xpos[3:4]), line=2.25, cex=0.6)

mtext("bidirectional,\nother", side=1, at=mean(xpos[5:6]), line=2.25, cex=0.6)
mtext("enhancer\noverlap", side=1, at=mean(xpos[7:8]), line=2.25, cex=0.6)

mtext("intergenic", side=1, at=mean(xpos[9]), line=1.25, cex=0.6)

mtext("promoters", side=3, line=0.15, cex=0.65)

legend("topright", legend=c("protein-coding", "lncRNA"), fill=c("indianred", "steelblue"), bty="n", inset=c(-0.01,-0.01), cex=0.95, horiz=F, xpd=NA)

mtext("B", side=3, at=-0.8, font=2, cex=0.9, line=0.75)

###################################################################################


## phastcons for splice sites

ylim=c(0,1)
xlim=c(0.5,3.5)

par(mar=c(5.1,3.1,2.2,0.5))

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

boxplot(splicesites.score[intersect(pc, spliced)], boxwex=0.7, col="indianred", at=1, add=T, axes=F, outline=F, notch=T)
boxplot(splicesites.score[intersect(lnc, spliced)], boxwex=0.7, col="steelblue", at=2, add=T, axes=F, outline=F, notch=T)

boxplot(intergenic.score, boxwex=0.7, col="gray60", at=3, add=T, axes=F, outline=F, notch=T)

n=c(length(intersect(pc, spliced)), length(intersect(lnc, spliced)))

mtext(n, side=1, cex=0.6, line=3.9, at=c(1,2))


axis(side=1, at=c(1:3), labels=rep("",3))
axis(side=2, cex.axis=0.85, mgp=c(3,0.5,0))

mtext("sequence conservation", side=2, line=2, cex=0.65)
mtext(c("protein\ncoding", "lncRNA", "intergenic"), side=1, at=1:3, line=0.75, cex=0.6, las=2)

mtext("splice sites", side=3, line=0.15, cex=0.65)

mtext("C", side=3, at=-0.75, font=2, cex=0.9, line=0.75)

###################################################################################

dev.off()

###################################################################################
