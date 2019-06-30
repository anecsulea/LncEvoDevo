##############################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE

  library(ape)

}

##########################################################################

if(load==TRUE){

  print("loading data")

  load("RData/data.expression.ortho.RData")
  load("RData/data.expression.conservation.spearman.RData")

  load=FALSE
}

##############################################################################

if(prepare==TRUE){

  ## example between-within species correlation

  lnc.mbr1=exportho.mr[which(exportho.mr$GeneType=="lncRNA"),"Mouse_Brain_EarlyEmbryo1"]
  lnc.mbr2=exportho.mr[which(exportho.mr$GeneType=="lncRNA"),"Mouse_Brain_EarlyEmbryo2"]
  lnc.mbr=(lnc.mbr1+lnc.mbr2)/2
  
  lnc.rbr1=exportho.mr[which(exportho.mr$GeneType=="lncRNA"),"Rat_Brain_EarlyEmbryo1"]
  lnc.rbr2=exportho.mr[which(exportho.mr$GeneType=="lncRNA"),"Rat_Brain_EarlyEmbryo2"]
  lnc.rbr=(lnc.rbr1+lnc.rbr2)/2

  
  pc.mbr1=exportho.mr[which(exportho.mr$GeneType=="protein_coding"),"Mouse_Brain_EarlyEmbryo1"]
  pc.mbr2=exportho.mr[which(exportho.mr$GeneType=="protein_coding"),"Mouse_Brain_EarlyEmbryo2"]
  pc.mbr=(pc.mbr1+pc.mbr2)/2

  pc.rbr1=exportho.mr[which(exportho.mr$GeneType=="protein_coding"),"Rat_Brain_EarlyEmbryo1"]
  pc.rbr2=exportho.mr[which(exportho.mr$GeneType=="protein_coding"),"Rat_Brain_EarlyEmbryo2"]
  pc.rbr=(pc.rbr1+pc.rbr2)/2
  
  prepare=FALSE
 }

###################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "Figure7.pdf", sep=""), width=5.25, height=5)

###################################################################################

m=matrix(rep(NA,10*8), nrow=10)

for(i in 1:4){
  m[i,]=c(rep(1:4, each=2))
}

for(i in 5:10){
  m[i,]=c(rep(5,4), rep(6,4))
}

layout(m)

###################################################################################

## example 

par(mar=c(2.5, 2.25, 2.1, 0))

## pc within species

smoothScatter(log2(pc.mbr1+1), log2(pc.mbr2+1), xlab="", ylab="", axes=F)
axis(side=2, cex.axis=0.8, mgp=c(3, 0.5, 0))
axis(side=1, cex.axis=0.8, mgp=c(3, 0.35, 0))
axis(side=1, cex.axis=0.8, at=12, mgp=c(3,0.35,0))
     
mtext("log2-transformed TPM", side=1, line=1.5, cex=0.65, at=32)
mtext("log2-transformed TPM", side=2, line=1.25, cex=0.65)

rho=cor(log2(pc.mbr1+1), log2(pc.mbr2+1), method="spearman", use="complete.obs")
mtext(paste("rho=", round(rho, digits=2), sep=""), side=3, cex=0.65)
text("within\nspecies", x=0, y=12, adj=0, cex=0.85)


mtext("A", side=3, at=-3.5, line=0.75, font=2, cex=0.95)

par(mar=c(2.5, 1.5, 2.1, 0.75))

## pc between species

smoothScatter(log2(pc.mbr+1), log2(pc.rbr+1), xlab="", ylab="", axes=F)

axis(side=2, cex.axis=0.8, mgp=c(3, 0.5, 0))
axis(side=1, cex.axis=0.8, mgp=c(3, 0.35, 0))
axis(side=1, cex.axis=0.8, at=12, mgp=c(3,0.35,0))

## mtext("mouse average", side=1, line=1, cex=0.6)
## mtext("rat average", side=2, line=1, cex=0.6)

rho=cor(log2(pc.mbr+1), log2(pc.rbr+1), method="spearman", use="complete.obs")
mtext(paste("rho=", round(rho, digits=2), sep=""), side=3, cex=0.65)
text("between\nspecies", x=0, y=11.5, adj=0, cex=0.85)

mtext("protein-coding genes", side=3, line=1, at=-1, cex=0.65)


## lnc within species
par(mar=c(2.5, 2.25, 2.1, 0))

smoothScatter(log2(lnc.mbr1+1), log2(lnc.mbr2+1), xlab="", ylab="", axes=F)
axis(side=2, cex.axis=0.8, mgp=c(3, 0.5, 0))
axis(side=1, cex.axis=0.8, mgp=c(3, 0.35, 0))

## mtext("mouse replicate 1", side=1, line=1, cex=0.6)
## mtext("mouse replicate 2", side=2, line=1, cex=0.6)

rho=cor(log2(lnc.mbr1+1), log2(lnc.mbr2+1), method="spearman", use="complete.obs")
mtext(paste("rho = ", round(rho, digits=2), sep=""), side=3, cex=0.65)
text("within\nspecies", x=0, y=6, adj=0, cex=0.85)

## lnc between species

par(mar=c(2.5, 1.5, 2.1, 0.75))
smoothScatter(log2(lnc.mbr+1), log2(lnc.rbr+1), xlab="", ylab="", axes=F)
axis(side=2, cex.axis=0.8, mgp=c(3, 0.5, 0))
axis(side=1, cex.axis=0.8, mgp=c(3, 0.35, 0))

## mtext("mouse average", side=1, line=1, cex=0.6)
## mtext("rat average", side=2, line=1, cex=0.6)

rho=cor(log2(lnc.mbr+1), log2(lnc.rbr+1), method="spearman", use="complete.obs")
mtext(paste("rho=", round(rho, digits=2), sep=""), side=3, cex=0.65)
text("between\nspecies", x=0, y=5.6, adj=0, cex=0.85)

mtext("lncRNAs", side=3, line=1, at=-1, cex=0.65)

###################################################################################


## plot expression conservation, pc genes

spacer=6
basicxpos=1:5
allxpos=c(basicxpos, basicxpos+spacer, basicxpos+2*spacer, basicxpos+3*spacer)
xpos=c(basicxpos, basicxpos+spacer, basicxpos+2*spacer, basicxpos[-1]+3*spacer)

this.tiss=unlist(lapply(names(corpc.between.tpm), function(x) unlist(strsplit(x, split="_"))[1]))
this.age=unlist(lapply(names(corpc.between.tpm), function(x) unlist(strsplit(x, split="_"))[2]))

## CI and yaxis

low.pc.tpm=apply(bootstrap.pc.between.tpm/bootstrap.pc.within.tpm, 2, min)
high.pc.tpm=apply(bootstrap.pc.between.tpm/bootstrap.pc.within.tpm, 2, max)

ylim=range(c(low.pc.tpm, high.pc.tpm))

## actual plot

par(mar=c(3.65, 3.15, 2.75, 0.65))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(min(xpos)-0.5, max(xpos)+0.5), ylim=ylim)

points(xpos, corpc.between.tpm/corpc.within.tpm, pch=pch.age[this.age], col=col.tissues[this.tiss], bg=col.tissues[this.tiss])
box()

meantiss=c(mean(xpos[1:5]), mean(xpos[6:10]),   mean(xpos[11:15]), mean(xpos[16:19]))

## axis(side=1, cex.axis=0.8, mgp=c(3,0.6, 0), at=meantiss, labels=rep("", 4))

axis(side=1, at=allxpos, labels=rep("", 20))
mtext(rep(1:5, 4), side=1, at=allxpos, line=0.5, cex=0.55)
mtext("dev. stage", side=1, line=0.5, at=-2.5, cex=0.55)


axis(side=2, cex.axis=0.85, mgp=c(3,0.65, 0))

mtext(shortname.tiss, side=1, at=meantiss, line=1.5, cex=0.65)
tinyx=(xpos[6]-xpos[5])/2

abline(v=c(mean(xpos[5:6]), mean(xpos[10:11]), xpos[15]+tinyx), lty=2)

segments(xpos, low.pc.tpm, xpos, high.pc.tpm, col=col.tissues[this.tiss], lwd=1.25)

mtext("expression conservation", side=2, line=1.75, cex=0.65)
mtext("protein-coding genes", side=3, line=0.25, cex=0.65)
mtext("B", side=3, at=-4, line=0.5, font=2, cex=0.95)

#############################################################################

## plot expression conservation, lncRNA genes

spacer=6
basicxpos=1:5
allxpos=c(basicxpos, basicxpos+spacer, basicxpos+2*spacer, basicxpos+3*spacer)
xpos=c(basicxpos, basicxpos+spacer, basicxpos+2*spacer, basicxpos[-1]+3*spacer)

this.tiss=unlist(lapply(names(corlnc.between.tpm), function(x) unlist(strsplit(x, split="_"))[1]))
this.age=unlist(lapply(names(corlnc.between.tpm), function(x) unlist(strsplit(x, split="_"))[2]))

## CI and yaxis

low.lnc.tpm=apply(bootstrap.lnc.between.tpm/bootstrap.lnc.within.tpm, 2, min)
high.lnc.tpm=apply(bootstrap.lnc.between.tpm/bootstrap.lnc.within.tpm, 2, max)

ylim=range(c(low.lnc.tpm, high.lnc.tpm))

## actual plot

par(mar=c(3.65, 3.15, 2.75, 0.65))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(min(xpos)-0.5, max(xpos)+0.5), ylim=ylim)

points(xpos, corlnc.between.tpm/corlnc.within.tpm, pch=pch.age[this.age], col=col.tissues[this.tiss], bg=col.tissues[this.tiss])
box()

meantiss=c(mean(xpos[1:5]), mean(xpos[6:10]),   mean(xpos[11:15]), mean(xpos[16:19]))

## axis(side=1, cex.axis=0.8, mgp=c(3,0.6, 0), at=meantiss, labels=rep("", 4))
axis(side=1, at=allxpos, labels=rep("", 20))
mtext(rep(1:5, 4), side=1, at=allxpos, line=0.5, cex=0.55)
mtext("dev. stage", side=1, line=0.5, at=-2.5, cex=0.55)

axis(side=2, cex.axis=0.85, mgp=c(3,0.65, 0))

mtext(shortname.tiss, side=1, at=meantiss, line=1.5, cex=0.65)
tinyx=(xpos[6]-xpos[5])/2

abline(v=c(mean(xpos[5:6]), mean(xpos[10:11]), xpos[15]+tinyx), lty=2)

segments(xpos, low.lnc.tpm, xpos, high.lnc.tpm, col=col.tissues[this.tiss], lwd=1.25)

mtext("expression conservation", side=2, line=1.75, cex=0.65)
mtext("lncRNAs", side=3, line=0.25, cex=0.65)
mtext("C", side=3, at=-4, line=0.5, font=2, cex=0.95)

###################################################################################

dev.off()

###################################################################################

