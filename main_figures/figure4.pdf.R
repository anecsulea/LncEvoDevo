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

pdf(file=paste(pathFigures, "Figure4.pdf", sep=""), width=6.85, height=5)

###################################################################################

## layout

m=matrix(rep(NA,2*30), nrow=2)

for(i in c(1)){
  m[i,]=rep(c(1,3,5), each=10)
}

for(i in c(2)){
  m[i,]=rep(c(2,4,6), each=10)
}

layout(m)

###################################################################################


## plot by types & max expression

labels=c("A", "B", "C")
names(labels)=c("exons", "promoters", "splicesites")
text=c("exons", "promoters", "splice sites")
names(text)=names(labels)

samples=kronecker(tissue.order, age.order, paste, sep="_")
xpos=1:20
xpos=setdiff(xpos, 16)

this.tissues=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))

for(type in c("exons", "promoters", "splicesites")){

  for(geneclass in c("pc", "lnc")){

    ylim=range(c(as.numeric(allvalues[[type]][[paste("ci.low.", geneclass, sep="")]]), as.numeric(allvalues[[type]][[paste("ci.high.", geneclass, sep="")]])))

    diffy=diff(ylim)

    ylim=ylim+c(-diffy/20, diffy/20)

    if(geneclass=="lnc"){
      if(ylim[1]>0.07){
        ylim[1]=0.07
      }
    }

    yax=pretty(ylim)

    ylim=c(min(c(ylim, yax)), max(c(ylim, yax)))
    
    ## yax=yax[which(yax>=ylim[1] & yax<=ylim[2])]
    
    xlim=c(0, 22)

    this.medians=as.numeric(allvalues[[type]][[paste("medians.", geneclass, sep="")]])
    this.low=as.numeric(allvalues[[type]][[paste("ci.low.", geneclass, sep="")]])
    this.high=as.numeric(allvalues[[type]][[paste("ci.high.", geneclass, sep="")]])
  
    if(geneclass=="pc"){
      par(mar=c(2.5, 3.1, 2.5, 1.1))
    } else{
      par(mar=c(3.5, 3.1, 1.5, 1.1))
    }

    plot(xpos, this.medians, type="n", xlim=xlim, ylim=ylim, axes=F, xaxs="i", yaxs="i", xlab="", ylab="")
        
    segments(xpos, ylim[1],xpos, ylim[2],lty=2, col="gray80")
    segments(xpos,this.low, xpos, this.high, col=col.tissues[this.tissues])
    points(xpos, this.medians, pch=20, col=col.tissues[this.tissues])
    
    ## intergenic
    segments(21, ci.low.ig, 21, ci.high.ig, col="gray60")
    segments(21, ylim[1], 21, ylim[2], lty=2, col="gray80")
    points(21, median.ig, pch=20, col="gray60")
    
    if(geneclass=="lnc"){
      axis(side=1, at=1:21, labels=rep("",21), mgp=c(3,0.75,0), cex.axis=0.85)
      mtext(rep(1:5, 4), side=1, at=1:20, cex=0.6, line=0.5)
      mtext("ig", side=1, at=21, cex=0.6, line=0.5, las=1)
      
      mtext("sequence conservation", side=2, line=2, cex=0.65)
      mtext("developmental stage", side=1, line=1.5, cex=0.65, at=mean(1:20))
      
    } else{
      mtext(text[type], side=3, cex=0.7)
      mtext(labels[type], side=3, at=-4.2, line=1, font=2, cex=0.9)
      mtext("sequence conservation", side=2, line=2, cex=0.65)
    }
    
        
    axis(side=2,  mgp=c(3,0.5,0), cex.axis=0.85, at=yax)

    if(type=="splicesites"){
      if(geneclass=="lnc"){
        mtext("lncRNAs", side=4, line=-0.35, cex=0.65, font=1)
      }
      
      if(geneclass=="pc"){
        mtext("protein-coding genes", side=4, line=-0.35, cex=0.65, font=1)
      }
    }

    if(type=="exons" & geneclass=="lnc"){

      legend("bottomleft", legend=c(shortname.tiss), fill=c(col.tissues), border=c(col.tissues), bty="n", inset=c(0.01,0.02), cex=0.95, horiz=F, xpd=NA, bg="white")

    }
  }
}


###################################################################################

dev.off()

###################################################################################
