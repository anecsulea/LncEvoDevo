##################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE
}

####################################################################

if(load==TRUE){

  print("loading data")

  load("RData/data.expression.divergence.RData")

  load("RData/data.aln.stats.RData")
  alnstat=alnstat.mr

  load("RData/data.expression.ortho.RData")
  rownames(exportho.mr)=exportho.mr$ID
  exportho.mr=exportho.mr[,which(!(colnames(exportho.mr)%in%c("ID", "GeneType")))]

  rownames(avgexp.mr)=avgexp.mr$ID
  avgexp.mr=avgexp.mr[,which(!(colnames(avgexp.mr)%in%c("ID", "GeneType")))]

  load("RData/data.annotations.Mouse.RData")
  allinfo.mouse=allinfo
  
  load=FALSE
}

######################################################################

if(prepare==TRUE){
 
  expdiv.pc=expdiv[which(expdiv$GeneType=="protein_coding"), "ExpressionDivergence"]
  names(expdiv.pc)=expdiv$ID[which(expdiv$GeneType=="protein_coding")]

  meantpm.pc=expdiv$MeanTPM[which(expdiv$GeneType=="protein_coding")]
  names(meantpm.pc)=names(expdiv.pc)
  
  expdiv.lnc=expdiv[which(expdiv$GeneType=="lncRNA"), "ExpressionDivergence"]
  names(expdiv.lnc)=expdiv$ID[which(expdiv$GeneType=="lncRNA")]
  
  meantpm.lnc=expdiv$MeanTPM[which(expdiv$GeneType=="lncRNA")]
  names(meantpm.lnc)=names(expdiv.lnc)
  
  alnstat$UngappedFraction.Mouse=alnstat$LengthUngapped/alnstat$ExonicLength.Mouse
  alnstat$UngappedFraction.Rat=alnstat$LengthUngapped/alnstat$ExonicLength.Rat
  
  alnstat$IdenticalFraction.Mouse=alnstat$LengthIdentical/alnstat$LengthUngapped
  alnstat$IdenticalFraction.Rat=alnstat$LengthIdentical/alnstat$LengthUngapped
  
  alnstat$MaxUngappedFraction=apply(alnstat[,c("UngappedFraction.Mouse", "UngappedFraction.Rat")],1,max)
  alnstat$MaxIdenticalFraction=apply(alnstat[,c("IdenticalFraction.Mouse", "IdenticalFraction.Rat")],1,max)

  alnstat$MinUngappedFraction=apply(alnstat[,c("UngappedFraction.Mouse", "UngappedFraction.Rat")],1,min)
  alnstat$MinIdenticalFraction=apply(alnstat[,c("IdenticalFraction.Mouse", "IdenticalFraction.Rat")],1,min)
    
  prepare=FALSE
}


######################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "Figure9.pdf", sep=""), width=6.85, height=5.5)

#############################################################################

## layout

m=matrix(rep(NA, 20*8), nrow=20)

for(i in 1:10){
  m[i,]=c(rep(1, 4), rep(2, 2), rep(3, 2))
}

for(i in 11:15){
  m[i,]=c(rep(4, 2), rep(5, 2), rep(6, 4))
}

for(i in 16:20){
  m[i,]=c(rep(4, 2), rep(5, 2), rep(7, 4))
}

layout(m)

######################################################################

## mean exp vs. expression divergence

genetype=c(rep("protein_coding", length(expdiv.pc)), rep("lncRNA", length(expdiv.lnc)))
names(genetype)=c(names(expdiv.pc), names(expdiv.lnc))

meanexp.allgenes=c(meantpm.pc, meantpm.lnc)
names(meanexp.allgenes)=names(genetype)

expdiv.allgenes=c(expdiv.pc, expdiv.lnc)
names(expdiv.allgenes)=names(genetype)


par(mar=c(3.5, 3.5, 2.1, 1.1))
smoothScatter(log2(meanexp.allgenes+1), expdiv.allgenes, xlab="", ylab="", axes=F)

mtext("average expression (log2 transformed TPM)", cex=0.65, side=1, line=1.75)
mtext("expression divergence, per gene", cex=0.65, side=2, line=2.25)

axis(side=2, cex.axis=0.85, mgp=c(3, 0.75, 0))
axis(side=1, cex.axis=0.85, mgp=c(3, 0.5, 0))

lm1=lm(expdiv.allgenes~log2(meanexp.allgenes+1))

r2=summary(lm1)$r.squared

abline(lm1, col="red")

text(paste("R2 =",round(r2,digits=2)), x=14, y=1.2, cex=0.95, xpd=NA)

mtext("A", side=3, line=0.8, at=-2.5, font=2, cex=0.9)

residuals=lm1$residuals

######################################################################

ylim=c(0, 1)

par(mar=c(3.1, 3.1, 2.1, 1.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(0.5, 2.5), ylim=ylim)

boxplot(expdiv.pc, at=1, add=T, col="indianred", axes=F, outline=F, boxwex=0.75)
boxplot(expdiv.lnc, at=2, add=T, col="steelblue", axes=F, outline=F, boxwex=0.75)

axis(side=2, mgp=c(3,0.75,0), cex.axis=0.75)
axis(side=1, mgp=c(3, 0.5, 0), at=c(1,2), labels=rep("",2), cex.axis=0.75)

mtext("expression divergence, per gene", side=2, cex=0.65, line=2.25)
mtext(c("pc", "lnc"), at=c(1,2), side=1, line=0.5, cex=0.65)

legend("topleft", legend=c("protein-coding", "lncRNA"), fill=c("indianred", "steelblue"), bty="n", inset=c(0.01, -0.01), cex=0.95, xpd=NA)

mtext("B", side=3, at=-0.22, font=2, cex=0.9, line=0.8)

######################################################################

ylim=c(-0.4, 0.6)

par(mar=c(3.1, 3.1, 2.1, 1.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(0.5, 2.5), ylim=ylim)

boxplot(residuals[which(genetype[names(residuals)]=="protein_coding")], at=1, add=T, col="indianred", axes=F, outline=F, boxwex=0.75)
boxplot(residuals[which(genetype[names(residuals)]=="lncRNA")], at=2, add=T, col="steelblue", axes=F, outline=F, boxwex=0.75)

axis(side=2, mgp=c(3,0.75,0), cex.axis=0.75)
axis(side=1, mgp=c(3, 0.5, 0), at=c(1,2), labels=rep("",2), cex.axis=0.75)

mtext("residual expression divergence", side=2, cex=0.65, line=2.25)
mtext(c("pc", "lnc"), at=c(1,2), side=1, line=0.5, cex=0.65)

mtext("C", side=3, at=-0.22, font=2, cex=0.9, line=0.8)

######################################################################

## expression divergence vs sequence divergence, protein-coding genes

par(mar=c(3.5, 3.5, 2.1, 0.1))
smoothScatter(100*alnstat[names(expdiv.pc),"MinUngappedFraction"], expdiv.pc, xlab="", ylab="", axes=F)

mtext("expression divergence, per gene", cex=0.65, side=2, line=2.25)


axis(side=2, cex.axis=0.85, mgp=c(3, 0.75, 0))
axis(side=1, cex.axis=0.85, mgp=c(3, 0.5, 0))

thisx=(100*alnstat[names(expdiv.pc),"MinUngappedFraction"])

lm2=lm(expdiv.pc~thisx)

r2=summary(lm2)$r.squared

abline(lm2, col="red")

text(paste("R2 =",round(r2,digits=2)), x=87, y=1.2, cex=0.95, xpd=NA)

mtext("protein-coding", cex=0.65, side=3, at=30, line=0.25)


mtext("D", side=3, line=0.8, at=-25, font=2, cex=0.9)

######################################################################

## expression divergence vs sequence divergence, lncRNAs

par(mar=c(3.5, 2, 2.1, 1.6))
smoothScatter(100*alnstat[names(expdiv.lnc),"MinUngappedFraction"], expdiv.lnc, xlab="", ylab="", axes=F)

axis(side=2, cex.axis=0.85, mgp=c(3, 0.75, 0))
axis(side=1, cex.axis=0.85, mgp=c(3, 0.5, 0))

thisx=(100*alnstat[names(expdiv.lnc),"MinUngappedFraction"])
lm3=lm(expdiv.lnc~thisx)

r2=summary(lm3)$r.squared

abline(lm3, col="red")

text(paste("R2 =",round(r2,digits=2)), x=77, y=1.2, cex=0.95, xpd=NA)

mtext("lncRNA", cex=0.65, side=3, at=12, line=0.25)

mtext("exonic sequence conservation (% ungapped sequence)", at=-0.1, side=1, cex=0.65, line=2)

######################################################################

## for each tissue/stage, expression divergence contribution

samples=kronecker(tissue.order, age.order, paste, sep="_")
xpos=1:20
names(xpos)=samples
xpos=setdiff(xpos, 16)
names(xpos)=setdiff(samples, "Testis_EarlyEmbryo")

xlim=c(0.5,20.5)
ylim=c(0,30)

## protein-coding

par(mar=c(0.9,2.75,2.1,1.1))
plot(1,type="n", xlab="", xlim=xlim, ylim=ylim, axes=F)

for(tiss in tissue.order){
  for(age in age.order){
    sample=paste(tiss, age, sep="_")

    if(paste("PCExpDiv",sample,sep=".")%in%colnames(expdiv)){
      this.expdiv.pc=expdiv[which(expdiv$GeneType=="protein_coding"), paste("PCExpDiv",sample,sep=".")]
           
      mean.pc=mean(this.expdiv.pc, na.rm=T)
     
      rect(xpos[sample]-0.2, 0, xpos[sample]+0.2, mean.pc, col=col.tissues[tiss], border=col.tissues[tiss])
    }
  }
}

axis(side=2, cex.axis=0.85)

legend("topleft", legend=shortname.tiss, border=col.tissues, fill=col.tissues, horiz=F, bty="n", xpd=NA, inset=c(0, -0.1), cex=1)
mtext("protein-coding", cex=0.65, side=3, at=11, line=0-.25)

mtext("E", side=3, at=-2.85, font=2, cex=0.9, line=0.8)

##################

## lncRNA

par(mar=c(3.0,2.75,0,1.1))
plot(1,type="n", xlab="", xlim=xlim, ylim=ylim, axes=F)

for(tiss in tissue.order){
  for(age in age.order){
    sample=paste(tiss, age, sep="_")

    if(paste("PCExpDiv",sample,sep=".")%in%colnames(expdiv)){
      this.expdiv.lnc=expdiv[which(expdiv$GeneType=="lncRNA"), paste("PCExpDiv",sample,sep=".")]
           
      mean.lnc=mean(this.expdiv.lnc, na.rm=T)
      rect(xpos[sample]-0.2, 0, xpos[sample]+0.2, mean.lnc, col=col.tissues[tiss], border=col.tissues[tiss])
  
    }
  }
}

axis(side=1, at=1:20, labels=rep("",20), mgp=c(3,0.5,0), cex.axis=0.7)
mtext(rep(1:5, 4), side=1, at=1:20, cex=0.5, line=0.5)

axis(side=2, cex.axis=0.85)

mtext("developmental stages", side=1, line=1.5, cex=0.65)
mtext("lncRNAs", cex=0.65, side=3, at=11, line=-0.25)

mtext("avg. part of expression divergence (%)", side=2, cex=0.65, at=35, line=2.25)


######################################################################

dev.off()

######################################################################
