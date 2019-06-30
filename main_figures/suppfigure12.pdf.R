##########################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE

}

##########################################################################

if(load==TRUE){

  load("RData/data.aln.stats.RData")
  alnstat=alnstat.mr

  load("RData/data.annotations.Mouse.RData")
  allinfo.mouse=allinfo
  pc.mouse=pc
  lnc.mouse=lnc

  load("RData/data.annotations.Rat.RData")
  allinfo.rat=allinfo
  pc.rat=pc
  lnc.rat=lnc

  

  load=FALSE
}
##########################################################################

if(prepare==TRUE){

  alnstat$FractionUngapped=alnstat$LengthUngapped/(apply(alnstat[,c("ExonicLength.Mouse", "ExonicLength.Rat")],1, max))

  alnstat$FractionIdentical=alnstat$LengthIdentical/alnstat$LengthUngapped

  alnstat$NbExons.Mouse=allinfo.mouse[alnstat$ID.Mouse, "NbExons"]
  alnstat$NbExons.Rat=allinfo.rat[alnstat$ID.Rat, "NbExons"]

  alnstat$DiffExons=abs(alnstat$NbExons.Mouse-alnstat$NbExons.Rat)/apply(alnstat[,c("NbExons.Mouse", "NbExons.Rat")],1, max)

  
  ## select alnstat
  alnstat.lnc=alnstat[which(alnstat$ID.Mouse%in%lnc.mouse & alnstat$ID.Rat%in%lnc.rat),]
  alnstat.pc=alnstat[which(alnstat$ID.Mouse%in%pc.mouse & alnstat$ID.Rat%in%pc.rat),]
  

  prepare=FALSE

}

##########################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

############################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure12.pdf", sep=""), width=6.85, height=2.75)

m=matrix(c(1,2,3), nrow=1, byrow=T)

layout(m)

######################################################################

## fraction ungapped length, pc and lncRNAs

d.pc=density(100*alnstat.pc$FractionUngapped, bw=2)
d.lnc=density(100*alnstat.lnc$FractionUngapped, bw=2)

xlim=range(d.pc$x, d.lnc$x)
ylim=range(d.pc$y, d.lnc$y)
ylim[2]=ylim[2]+diff(ylim)/10

par(mar=c(3.5,3.1,2.1,1))

plot(1, type="n", xlab="", ylab="", ylim=ylim, xlim=xlim, axes=F)

axis(side=1, cex.axis=0.85, mgp=c(3,0.5,0))
axis(side=2, cex.axis=0.85, mgp=c(3,0.5,0))

mtext("% aligned exonic sequence", side=1, line=1.75, cex=0.7)
mtext("density", side=2, line=1.75, cex=0.7)

lines(d.pc, col="indianred", lty=1)
lines(d.lnc, col="steelblue", lty=1)

mtext("A", side=3, at=-30, line=1, font=2, cex=0.9)

######################################################################

## fraction identical length, pc and lncRNAs

d.pc=density(100*alnstat.pc$FractionIdentical, bw=1)
d.lnc=density(100*alnstat.lnc$FractionIdentical, bw=1)

xlim=range(d.pc$x, d.lnc$x)
ylim=range(d.pc$y, d.lnc$y)
ylim[2]=ylim[2]+diff(ylim)/10

par(mar=c(3.5,3.1,2.1,1))

plot(1, type="n", xlab="", ylab="", ylim=ylim, xlim=xlim, axes=F)

axis(side=1, cex.axis=0.85, mgp=c(3,0.5,0))
axis(side=2, cex.axis=0.85, mgp=c(3,0.5,0))

mtext("% identical aligned sequence", side=1, line=1.5, cex=0.65)
mtext("density", side=2, line=1.75, cex=0.65)

lines(d.pc, col="indianred", lty=1)
lines(d.lnc, col="steelblue", lty=1)

mtext("B", side=3, at=48, line=1, font=2, cex=0.9)

######################################################################

## diff nbexons , pc and lncRNAs

d.pc=density(100*alnstat.pc$DiffExons, bw=1)
d.lnc=density(100*alnstat.lnc$DiffExons, bw=1)

xlim=range(d.pc$x, d.lnc$x)
ylim=range(d.pc$y, d.lnc$y)
ylim[2]=ylim[2]+diff(ylim)/10

par(mar=c(3.5,3.1,2.1,1))

plot(1, type="n", xlab="", ylab="", ylim=ylim, xlim=xlim, axes=F)

axis(side=1, cex.axis=0.85, mgp=c(3,0.5,0))
axis(side=2, cex.axis=0.85, mgp=c(3,0.5,0))

mtext("relative difference nb. exons", side=1, line=1.5, cex=0.65)
mtext("density", side=2, line=1.75, cex=0.65)

lines(d.pc, col="indianred", lty=1)
lines(d.lnc, col="steelblue", lty=1)

legend("topright", c("protein-coding", "lncRNA"), lty=1, col=c("indianred", "steelblue"), bty="n")


mtext("C", side=3, at=-25, line=1, font=2, cex=0.9)

######################################################################


######################################################################

dev.off()

######################################################################
