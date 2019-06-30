###################################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE
}

###################################################################################

if(load==TRUE){
  for(sp in c("Mouse", "Rat", "Chicken")){
    load(paste("RData/data.annotations.", sp, ".RData", sep=""))
    assign(paste("pc", tolower(sp), sep="."), pc)
    assign(paste("lnc", tolower(sp), sep="."), lnc)
    assign(paste("allinfo", tolower(sp), sep="."), allinfo)
    
    load(paste("RData/data.expression.",sp,".RData", sep=""))
    assign(paste("downsampled", tolower(sp), sep="."), downsampled)
    
    rm(list=c("pc","lnc", "allinfo", "rawtpm", "normtpm", "kcounts", "readcounts", "downsampled"))
  }
  
  load=F
}

###################################################################################

if(prepare==TRUE){
  ## nb detected genes in each tissue and age

  minreads=10

  ## pc 
  nb.expressed.pc=list()
  
  nb.expressed.pc[["mouse"]]=apply(downsampled.mouse[pc.mouse,],2,function(x) length(which(x>=minreads)))
  names(nb.expressed.pc[["mouse"]])=colnames(downsampled.mouse)
  
  nb.expressed.pc[["rat"]]=apply(downsampled.rat[pc.rat,],2,function(x) length(which(x>=minreads)))
  names(nb.expressed.pc[["rat"]])=colnames(downsampled.rat)
  
  nb.expressed.pc[["chicken"]]=apply(downsampled.chicken[pc.chicken,],2,function(x) length(which(x>=minreads)))
  names(nb.expressed.pc[["chicken"]])=colnames(downsampled.chicken)

  ## lncRNAs 
  nb.expressed.lnc=list()
  
  nb.expressed.lnc[["mouse"]]=apply(downsampled.mouse[lnc.mouse,],2,function(x) length(which(x>=minreads)))
  names(nb.expressed.lnc[["mouse"]])=colnames(downsampled.mouse)
  
  nb.expressed.lnc[["rat"]]=apply(downsampled.rat[lnc.rat,],2,function(x) length(which(x>=minreads)))
  names(nb.expressed.lnc[["rat"]])=colnames(downsampled.rat)
  
  nb.expressed.lnc[["chicken"]]=apply(downsampled.chicken[lnc.chicken,],2,function(x) length(which(x>=minreads)))
  names(nb.expressed.lnc[["chicken"]])=colnames(downsampled.chicken)
  
  prepare=F

}

###################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "Figure2.pdf", sep=""), width=5.5, height=2.75)

###################################################################################

## layout

m=matrix(rep(NA,5*12), nrow=5)

for(i in 1:5){
  m[i,]=c(rep(1,6), rep(2,6))
}

layout(m)

###################################################################################

## number of detected protein-coding genes

par(tck=NA)

yax=pretty(c(min(unlist(lapply(nb.expressed.pc, min))), max(unlist(lapply(nb.expressed.pc, max)))))
ylim=range(yax)
ylim[1]=min(11500, ylim[1])

xlim=c(0.6,5.4) ## 5 stages

par(mar=c(3.5, 3.5, 2, 1))

plot(1, type="n", axes=F, xlab="", ylab="", xlim=xlim, ylim=ylim)


xpos=1:5
names(xpos)=c("EarlyEmbryo", "LateEmbryo", "Newborn", "Adult", "Aged")

smallx=c(-0.25,0,0.25)
names(smallx)=c("mouse", "rat", "chicken")

pch.sp=c(21, 23, 24)
names(pch.sp)=c("mouse", "rat", "chicken")

for(tissue in tissue.order){
  this.col=col.tissues[tissue]

  for(age in age.order){
    for(sp in names(pch.sp)){
      this.x=xpos[age]+smallx[sp]
      this.pch=pch.sp[sp]
      
      points(this.x, nb.expressed.pc[[sp]][paste(tissue, age, sep="_")], pch=this.pch, col="black", bg=this.col, cex=1.35, lwd=1.25)
    }
  }
}

box()
axis(side=1, at=xpos, labels=rep("",5), mgp=c(3, 0.75, 0))


mtext(1:5,  at=xpos, line=0.45, cex=0.55, side=1)

axis(side=2, at=yax, labels=rep("",length(yax)), mgp=c(3,0.75,0))
mtext(paste(yax/1000, "K", sep=""), at=yax, side=2, cex=0.5, line=0.6)

mtext("# detected genes", side=2, line=2, cex=0.65, xpd=NA)
mtext("developmental stages", side=1, line=1.8, cex=0.65)

mtext("protein-coding genes", side=3, font=2, cex=0.7, line=0.15)

mtext("A", side=3, at=-0.3, font=2, cex=0.95, line=0.8)


#######################################################################################

## number of detected lncRNAs

ylim=c(min(unlist(lapply(nb.expressed.lnc, min))), max(unlist(lapply(nb.expressed.lnc, max))))
yax=pretty(ylim)

par(mar=c(3.5, 3.5, 2, 1))

plot(1, type="n", axes=F, xlab="", ylab="", xlim=xlim, ylim=range(c(yax, ylim)))

for(tissue in tissue.order){
  this.col=col.tissues[tissue]

  for(age in age.order){
    for(sp in names(pch.sp)){
      this.x=xpos[age]+smallx[sp]
      this.pch=pch.sp[sp]
      
      points(this.x, nb.expressed.lnc[[sp]][paste(tissue, age, sep="_")], pch=this.pch, col="black", bg=this.col, cex=1.35, lwd=1.25)
    }
  }
}

box()
axis(side=1, at=xpos, labels=rep("",5), mgp=c(3, 0.75, 0))

mtext(1:5,  at=xpos, line=0.5, cex=0.55, side=1)

axis(side=2, at=yax, labels=rep("",length(yax)), mgp=c(3,0.75,0))
mtext(paste(yax/1000, "K", sep=""), at=yax, side=2, cex=0.5, line=0.6)

mtext("# detected genes", side=2, line=2, cex=0.65,xpd=NA)
mtext("developmental stages", side=1, line=1.8, cex=0.65)

mtext("lncRNAs", side=3, font=2, cex=0.7, line=0.15)
mtext("B", side=3, at=-0.3, font=2, cex=0.95, line=0.8)

legend("topleft", legend=c("mouse", "rat", "chicken"), pch=pch.allsp, inset=c(0.01,-0.01), bty="n", cex=0.95, horiz=F, xpd=NA)

legend("topleft", legend=shortname.tiss, fill=col.tissues, border="black", bty="n", inset=c(-0.01,0.2), cex=0.95, horiz=F, xpd=NA)


#####################################################################################


dev.off()

###################################################################################
