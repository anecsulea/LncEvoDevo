################################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE
  
  maxFDR=0.01
  
  library(ape)
}

################################################################################

if(load==TRUE){

  for(sp in c("Mouse", "Rat", "Chicken")){
    load(paste("RData/data.annotations.", sp, ".RData", sep=""))
    assign(paste("pc", tolower(sp), sep="."), pc)
    assign(paste("lnc", tolower(sp), sep="."), lnc)
    assign(paste("allinfo", tolower(sp), sep="."), allinfo)
    
    load(paste("RData/data.tpm.stats.",sp,".RData", sep=""))
    assign(paste("avgtpm", tolower(sp), sep="."), stats)

    if(sp!="Chicken"){
      load(paste("RData/data.diffexp.",sp,".RData", sep=""))
      assign(paste("de.allreads.", tolower(sp), sep=""), de.global.allreads)
      assign(paste("de.resampled.", tolower(sp), sep=""), de.global.resampled)
    }
  }

  load=FALSE
}

################################################################################

if(prepare==TRUE){

  table.maxstage=list()
  ratio.signif=list()

  for(sp in c("Mouse", "Rat", "Chicken")){
    avgtpm=get(paste("avgtpm", tolower(sp), sep="."))

    this.info=get(paste("allinfo", tolower(sp), sep="."))
       
    this.pc=get(paste("pc", tolower(sp), sep="."))
    this.lnc=get(paste("lnc", tolower(sp), sep="."))

    ## we take all genes
    
    maxsample=avgtpm[,"MaxSample"]
    maxtissue=unlist(lapply(maxsample, function(x) unlist(strsplit(x, split="_"))[1]))
    maxage=unlist(lapply(maxsample, function(x) unlist(strsplit(x, split="_"))[2]))

    names(maxsample)=rownames(avgtpm)
    names(maxtissue)=rownames(avgtpm)
    names(maxage)=rownames(avgtpm)

    t.pc=as.numeric(table(as.factor(maxtissue[this.pc])))
    names(t.pc)=levels(as.factor(maxtissue))
    t.pc=100*t.pc/sum(t.pc)

    t.lnc=as.numeric(table(as.factor(maxtissue[this.lnc])))
    names(t.lnc)=levels(as.factor(maxtissue))
    t.lnc=100*t.lnc/sum(t.lnc)

    a.pc=as.numeric(table(factor(maxage[this.pc], levels=age.order)))
    names(a.pc)=age.order
    a.pc=100*a.pc/sum(a.pc)

    a.lnc=as.numeric(table(factor(maxage[this.lnc], levels=age.order)))
    names(a.lnc)=age.order
    a.lnc=100*a.lnc/sum(a.lnc)

      
    assign(paste("maxage", tolower(sp), sep="."), maxage)
    assign(paste("maxtissue", tolower(sp), sep="."), maxtissue)
   
    assign(paste("t.pc", tolower(sp), sep="."), t.pc)
    assign(paste("t.lnc", tolower(sp), sep="."), t.lnc)
    
    assign(paste("a.pc", tolower(sp), sep="."), a.pc)
    assign(paste("a.lnc", tolower(sp), sep="."), a.lnc)

    
    if(sp%in%c("Mouse", "Rat")){
      ratio.signif[[sp]]=list()
      
      for(type in c("allreads", "resampled")){
        ratio.signif[[sp]][[type]]=list()
        
        this.de=get(paste("de", type, tolower(sp), sep="."))
        
        this.signif.pc=unlist(lapply(tissue.order, function(x) length(which(this.de[,"GeneID"]%in%this.pc & this.de[,paste("FDR", x, sep=".")]<maxFDR))))
        this.signif.lnc=unlist(lapply(tissue.order, function(x) length(which(this.de[,"GeneID"]%in%this.lnc & this.de[,paste("FDR", x, sep=".")]<maxFDR))))
        
        this.all.pc=unlist(lapply(tissue.order, function(x) length(which(this.de[,"GeneID"]%in%this.pc & !is.na(this.de[,paste("FDR", x, sep=".")])))))
        this.all.lnc=unlist(lapply(tissue.order, function(x) length(which(this.de[,"GeneID"]%in%this.lnc & !is.na(this.de[,paste("FDR", x, sep=".")])))))
        
        ratio.signif[[sp]][[type]][["pc"]]=this.signif.pc/this.all.pc
        ratio.signif[[sp]][[type]][["lnc"]]=this.signif.lnc/this.all.lnc

      }

      ## table max stage in each species

      table.maxstage[[sp]]=list()

      this.de=get(paste("de.allreads", tolower(sp), sep="."))

      for(tiss in tissue.order){
        this.age.order=age.order

        if(tiss=="Testis"){
          this.age.order=setdiff(age.order, "EarlyEmbryo")
        }
        
        this.tpm=avgtpm[,paste("MeanTPM.",tiss, "_", this.age.order,sep="")]
        colnames(this.tpm)=this.age.order
        this.maxstage=apply(this.tpm,1, function(x) this.age.order[which.max(x)])
        names(this.maxstage)=rownames(this.tpm)
        
        this.signif=this.de$GeneID[which(this.de[,paste("FDR", tiss, sep=".")]<maxFDR)]
        this.signif.pc=intersect(this.pc, this.signif)
        this.signif.lnc=intersect(this.lnc, this.signif)

        this.table.pc=as.numeric(table(factor(this.maxstage[this.signif.pc], levels=age.order)))
        names(this.table.pc)=age.order

        this.table.pc=100*this.table.pc/sum(this.table.pc)

        this.table.lnc=as.numeric(table(factor(this.maxstage[this.signif.lnc], levels=age.order)))
        names(this.table.lnc)=age.order

        this.table.lnc=100*this.table.lnc/sum(this.table.lnc)

        table.maxstage[[sp]][[tiss]]=list("pc"=this.table.pc, "lnc"=this.table.lnc)
      }
    }
  }
  
  
  prepare=FALSE
}

################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in 
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "Figure3.pdf", sep=""), width=6.85, height=5.5)

m=matrix(rep(NA,10*10), nrow=10)

for(i in 1:5){
  m[i,]=c(rep(1, 3), rep(2, 3), rep(3, 4))
}

for(i in 6:10){
  m[i,]=c(rep(4,10))
}

layout(m)

################################################################################

## first plot: max tissue for pc and lnc, mouse and rat

par(mar=c(4.5, 3.1, 2.1, 0.5))

ylim=c(0, 101)
xlim=c(0.5, 7.5)
xpos=c(1, 2, 3.5, 4.5, 6, 7)
width=0.4

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

rect(xpos[1]-width, c(0, cumsum(t.pc.mouse)[1:3]), xpos[1]+width, cumsum(t.pc.mouse), col=col.tissues, border=col.tissues)
rect(xpos[2]-width, c(0, cumsum(t.lnc.mouse)[1:3]), xpos[2]+width, cumsum(t.lnc.mouse), col=col.tissues, border=col.tissues)

rect(xpos[3]-width, c(0, cumsum(t.pc.rat)[1:3]), xpos[3]+width, cumsum(t.pc.rat), col=col.tissues, border=col.tissues)
rect(xpos[4]-width, c(0, cumsum(t.lnc.rat)[1:3]), xpos[4]+width, cumsum(t.lnc.rat), col=col.tissues, border=col.tissues)

rect(xpos[5]-width, c(0, cumsum(t.pc.chicken)[1:2]), xpos[5]+width, cumsum(t.pc.chicken), col=col.tissues[1:3], border=col.tissues)
rect(xpos[6]-width, c(0, cumsum(t.lnc.chicken)[1:2]), xpos[6]+width, cumsum(t.lnc.chicken), col=col.tissues[1:3], border=col.tissues)

axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("",6))
mtext(rep("pc", 3), at=xpos[c(1,3,5)], side=1, line=0.52, cex=0.65)
mtext(rep("lnc", 3), at=xpos[c(2,4,6)], side=1, line=0.52, cex=0.65)
mtext(c("mouse", "rat", "chicken"), at=c(1.5, 4, 6.5), side=1, line=1.5, cex=0.65)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.9)

mtext("% genes w. max. expression in organ", side=2, cex=0.65, line=2.15)
mtext("A", side=3, at=-1.3, line=1, font=2, cex=0.95)


legend("topleft", legend=c("br", "kd", "lv", "ts"), fill=col.tissues, border=col.tissues, bty="n", inset=c(0.01,-0.1), cex=0.95, horiz=T, xpd=NA)

################################################################################

## second plot: max age for pc and lnc, mouse and rat

par(mar=c(4.5, 3.1, 2.1, 0.5))

ylim=c(0, 101)
xlim=c(0.5, 7.5)
xpos=c(1, 2, 3.5, 4.5, 6, 7)
width=0.4

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

rect(xpos[1]-width, c(0, cumsum(a.pc.mouse)[1:4]), xpos[1]+width, cumsum(a.pc.mouse), col=col.ages, border=col.ages)
rect(xpos[2]-width, c(0, cumsum(a.lnc.mouse)[1:4]), xpos[2]+width, cumsum(a.lnc.mouse), col=col.ages, border=col.ages)

rect(xpos[3]-width, c(0, cumsum(a.pc.rat)[1:4]), xpos[3]+width, cumsum(a.pc.rat), col=col.ages, border=col.ages)
rect(xpos[4]-width, c(0, cumsum(a.lnc.rat)[1:4]), xpos[4]+width, cumsum(a.lnc.rat), col=col.ages, border=col.ages)

rect(xpos[5]-width, c(0, cumsum(a.pc.chicken)[1]), xpos[5]+width, cumsum(a.pc.chicken[1:2]), col=col.ages[1:2], border=col.ages)
rect(xpos[6]-width, c(0, cumsum(a.lnc.chicken)[1]), xpos[6]+width, cumsum(a.lnc.chicken[1:2]), col=col.ages[1:2], border=col.ages)

axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("",6))
mtext(rep("pc", 3), at=xpos[c(1,3,5)], side=1, line=0.52, cex=0.65)
mtext(rep("lnc", 3), at=xpos[c(2,4,6)], side=1, line=0.52, cex=0.65)
mtext(c("mouse", "rat", "chicken"), at=c(1.5, 4, 6.5), side=1, line=1.5, cex=0.65)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.9)

mtext("% genes w. max. expression in stage", side=2, cex=0.65, line=2.15)
mtext("B", side=3, at=-1.3, line=1, font=2, cex=0.95)


legend("topleft", legend=c(1:5), fill=col.ages, border=col.ages, bty="n", inset=c(-0.05,-0.1), cex=0.95, horiz=T, xpd=NA)

################################################################################
## % of genes that are differentially expressed with time

xpos=1:4
names(xpos)=tissue.order

smallx=c(-0.2, 0.2)
names(smallx)=c("Mouse", "Rat")

bigx=c(0, 4.5)
names(bigx)=c("allreads", "resampled")

ylim=c(-5,105)
xlim=c(0, 10.5)

par(mar=c(3.75,3.5,1.1,0.1))

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i")

for(sp in c("Mouse", "Rat")){
  for(type in c("allreads", "resampled")){
    points(xpos+smallx[sp]+bigx[type], 100*ratio.signif[[sp]][[type]][["pc"]], pch=pch.allsp[tolower(sp)], col=col.tissues, bg=col.tissues)
    points(xpos+smallx[sp]+bigx[type], 100*ratio.signif[[sp]][[type]][["lnc"]], pch=pch.allsp[tolower(sp)], bg="white", col=col.tissues)
  }
}

abline(v=mean(c(xpos[4], xpos[1]+bigx["resampled"])), lty=2)

axis(side=1, at=c(xpos+bigx["allreads"], xpos+bigx["resampled"]), labels=rep("", 8), cex.axis=0.8, mgp=c(3, 0.75, 0))
axis(side=2, at=seq(from=0, to=100, by=20), cex.axis=0.9, mgp=c(3, 0.5, 0))

mtext("all reads", side=1, at=mean(xpos+bigx["allreads"]), line=1.5, cex=0.65)
mtext("downsampling", side=1, at=mean(xpos+bigx["resampled"]), line=1.5, cex=0.65)
mtext("% diff. exp. among stages", side=2, line=2, cex=0.65)

mtext(rep(c("br", "kd", "lv", "ts"), 2), at=c(xpos+bigx["allreads"], xpos+bigx["resampled"]), col=col.tissues, side=1, cex=0.65, line=0.45)

legend("bottomleft", legend=c("mouse", "rat"), pch=pch.allsp[c("mouse", "rat")], inset=-0.01, bty="n", cex=0.95)
legend("topright", legend=c("protein-coding", "lncRNA"), pch=22, col="black", pt.bg=c("black", "white"), inset=c(0.01, -0.05), bty="n", cex=0.95, xpd=NA)

mtext("C", side=3, line=0, at=-1.55, font=2, cex=0.95)

################################################################################

## maxstage for diffexp pc and lnc, for mouse only

this.table.maxstage=table.maxstage[["Mouse"]]

addspace=0.5
xpos=c(1:5, 6:10+addspace, 11:15+2*addspace, 16:20+3*addspace)
names(xpos)=kronecker(tissue.order, age.order, paste, sep="_")

ylim=c(0,60)
xlim=c(min(xpos)-0.5,max(xpos)+0.5)

width=0.2
space=0.04

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", xaxs="i")

for(tiss in tissue.order){
  for(age in age.order){
    this.x=xpos[paste(tiss,age, sep="_")]
    this.prop.pc=this.table.maxstage[[tiss]][["pc"]][age]
    this.prop.lnc=this.table.maxstage[[tiss]][["lnc"]][age]

    if(this.prop.pc>0 & this.prop.lnc>0){
      rect(this.x-width-space, 0, this.x-space, this.prop.pc, col=col.tissues[tiss], border=col.tissues[tiss])
      rect(this.x+space, 0, this.x+width+space, this.prop.lnc, col="white", border=col.tissues[tiss])
    }
  }
}

axis(side=2,  at=seq(from=0, to=60, by=20), cex.axis=0.9, mgp=c(3, 0.5, 0))
axis(side=1, at=xpos, labels=rep("",20), mgp=c(3,0.5,0))
mtext(rep(1:5, 4), at=xpos, line=0.5, cex=0.65, side=1)
mtext("% differentially expressed genes", side=2, line=2, cex=0.65)
mtext("developmental stage with maximum expression, mouse genes", side=1, line=2, cex=0.65)

abline(v=c(mean(xpos[5:6]), mean(xpos[10:11]), mean(xpos[15:16])), lty=2, col="gray40")

legend("topleft", legend=c("protein-coding", "lncRNA"), pch=22, col="black", pt.bg=c("black", "white"), inset=c(0.01, -0.05), bty="n", cex=0.95, xpd=NA)


mtext("D", side=3, line=0.5, at=-0.6, font=2, cex=0.95)

################################################################################

dev.off()

################################################################################




