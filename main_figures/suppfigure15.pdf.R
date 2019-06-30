##################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")

  source("plot.coverage.R")
  source("plot.annotations.R")
  
  load=TRUE
  prepare=TRUE
}

####################################################################

if(load==TRUE){
  load("RData/data.mouse.specific.candidate.RData")

  load("RData/data.annotations.Mouse.RData")
  allinfo.mouse=allinfo
  pc.mouse=pc
  rownames(allinfo.mouse)=allinfo.mouse$GeneID
  
  load("RData/data.annotations.Rat.RData")
  allinfo.rat=allinfo
  pc.rat=pc
  rownames(allinfo.rat)=allinfo.rat$GeneID
  
  load("RData/data.expression.divergence.RData")
  rownames(expdiv)=expdiv$ID
  
  load("RData/data.comparison.specific.ortho.RData")

  load=FALSE
}
####################################################################


if(prepare==TRUE){


  ## lm expression divergence

  expdiv=expdiv[which(expdiv$MeanTPM>0),]
  this.lm=lm(expdiv$ExpressionDivergence~log2(expdiv$MeanTPM+1))
  expdiv$ResidualExpressionDivergence=this.lm$residuals
  
  
  colnames(annot.mouse)=c("geneid", "exonid", "chr", "start", "end", "strand")

  mouse.genes=allinfo.mouse$GeneID[which(allinfo.mouse$Chr==mouse.chr & (((allinfo.mouse$Start>=mouse.start) & (allinfo.mouse$Start<=mouse.end)) | ((allinfo.mouse$End>=mouse.start) & (allinfo.mouse$End<=mouse.end))))]
  
  annot.mouse=annot.mouse[which(annot.mouse$geneid%in%mouse.genes),]
  annot.mouse$biotype=allinfo.mouse[annot.mouse$geneid, "EnsemblBiotype"]
  annot.mouse$genename=allinfo.mouse[annot.mouse$geneid, "GeneName"]

  colnames(annot.rat)=c("geneid", "exonid", "chr", "start", "end", "strand")
  
  rat.genes=allinfo.rat$GeneID[which(allinfo.rat$Chr==rat.chr & (((allinfo.rat$Start>=rat.start) & (allinfo.rat$Start<=rat.end)) | ((allinfo.rat$End>=rat.start) & (allinfo.rat$End<=rat.end))))]
  
  annot.rat=annot.rat[which(annot.rat$geneid%in%rat.genes),]
  
  annot.rat$biotype=allinfo.rat[annot.rat$geneid, "EnsemblBiotype"]
  annot.rat$genename=allinfo.rat[annot.rat$geneid, "GeneName"]

  ## neighbours

  neighbours.specific.mouse=intersect(unique(unlist(lapply(allinfo.mouse[mouse.specific, "BidirectionalPromoter"], function(x) unlist(strsplit(x, split=","))))),pc.mouse)

  neighbours.ortho.mouse=intersect(unique(unlist(lapply(allinfo.mouse[mouse.ortho, "BidirectionalPromoter"], function(x) unlist(strsplit(x, split=","))))),pc.mouse)

  common=intersect(neighbours.specific.mouse, neighbours.ortho.mouse)

  if(length(common)>0){
    neighbours.specific.mouse=setdiff(neighbours.specific.mouse, common)
    neighbours.ortho.mouse=setdiff(neighbours.ortho.mouse, common)
  }

  ortho.neighbours.specific.mouse=setdiff(unlist(lapply(neighbours.specific.mouse, function(x) grep(paste("^", x,"_",sep=""), expdiv$ID, value=T))), NA)
  ortho.neighbours.ortho.mouse=setdiff(unlist(lapply(neighbours.ortho.mouse, function(x) grep(paste("^", x,"_",sep=""), expdiv$ID, value=T))), NA)


  ## same for rat
  neighbours.specific.rat=intersect(unique(unlist(lapply(allinfo.rat[rat.specific, "BidirectionalPromoter"], function(x) unlist(strsplit(x, split=","))))),pc.rat)

  neighbours.ortho.rat=intersect(unique(unlist(lapply(allinfo.rat[rat.ortho, "BidirectionalPromoter"], function(x) unlist(strsplit(x, split=","))))),pc.rat)

  common=intersect(neighbours.specific.rat, neighbours.ortho.rat)

  if(length(common)>0){
    neighbours.specific.rat=setdiff(neighbours.specific.rat, common)
    neighbours.ortho.rat=setdiff(neighbours.ortho.rat, common)
  }

  ortho.neighbours.specific.rat=setdiff(unlist(lapply(neighbours.specific.rat, function(x) grep(paste("_", x,"$",sep=""), expdiv$ID, value=T))), NA)
  ortho.neighbours.ortho.rat=setdiff(unlist(lapply(neighbours.ortho.rat, function(x) grep(paste("_", x,"$",sep=""), expdiv$ID, value=T))), NA)

  ## combine both species

  ortho.neighbours.ortho=unique(c(ortho.neighbours.ortho.mouse,ortho.neighbours.ortho.rat))
  ortho.neighbours.specific=unique(c(ortho.neighbours.specific.mouse,ortho.neighbours.specific.rat))
  
  
  prepare=FALSE
}

####################################################################


## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure15.pdf", sep=""), width=4.49, height=6)

###################################################################################

## layout

m=matrix(rep(NA, 15*6), nrow=15)

for(i in c(1)){
  m[i,]=c(rep(1, 6)) ## mouse annotations
}

for(i in c(2:3)){
  m[i,]=c(rep(2, 6)) ## mouse kidney adult
}

for(i in c(4:5)){
  m[i,]=c(rep(3, 6)) ## mouse kidney aged
}


for(i in c(6)){
  m[i,]=c(rep(4, 6)) ## rat annotations
}

for(i in c(7:8)){
  m[i,]=c(rep(5, 6)) ## rat kidney adult
}

for(i in c(9:10)){
  m[i,]=c(rep(6, 6)) ## rat kidney aged
}

for(i in c(11:15)){
  m[i,]=c(rep(7,3), rep(8,3))
}

layout(m)
    
    
#####################################################################################

par(mar=c(0.1, 3.1, 0.5, 1.5)) ## for annotation plots


#####################################################################################
## mouse annotations
plot.annot.flat(annot=annot.mouse, start=mouse.start, end=mouse.end, strands=c("1"), biotypes=unique(annot.mouse$biotype), ylim=c(0.3,0.7), ypos=c(0.4,0.6), col.fwd="steelblue")

rect(proj.mouse.start, 0.45, proj.mouse.end, 0.55, border="red", col="white")

mtext("A", side=3, line=-1, at=mouse.start-(mouse.end-mouse.start)/15, font=2, cex=0.9)

#####################################################################################


par(mar=c(2.1, 3.1, 0.1, 1.5)) ## for coverage plots

for(tiss in c("Kidney_Adult", "Kidney_Aged")){
  plot.coverage(cov.mouse[[tiss]], xstart=mouse.start, xend=mouse.end, col=col.tissues[["Kidney"]],recompute.ylim=F, ylim=c(0,2000))

  if(tiss=="Kidney_Adult"){
    text="mouse kidney young adult"

    mtext("unique read coverage", at=-500, side=2, line=1.5, cex=0.65)

  } 
  if(tiss=="Kidney_Aged"){
    text="mouse kidney aged adult"
  }
   
  text(text, x=mouse.start+(mouse.end-mouse.start)/50, y=2000, cex=0.85, adj=c(0,1))
}

#####################################################################################

## rat annotations


par(mar=c(0.1, 3.1, 0.5, 1.5)) ## for annotation plots


plot.annot.flat(annot=annot.rat, start=rat.start, end=rat.end, strands=c("1"), biotypes=unique(annot.rat$biotype), ylim=c(0.35,0.65), ypos=c(0.5), col.fwd="steelblue")

rect(proj.rat.start, 0.45, proj.rat.end, 0.55, border="red", col="white")

#####################################################################################


par(mar=c(2.1, 3.1, 0.1, 1.5)) ## for coverage plots

for(tiss in c("Kidney_Adult", "Kidney_Aged")){
  plot.coverage(cov.rat[[tiss]], xstart=rat.start, xend=rat.end, col=col.tissues[["Kidney"]],recompute.ylim=F, ylim=c(0,2000))

  if(tiss=="Kidney_Adult"){
    text="rat kidney young adult"

    mtext("unique read coverage", at=-500, side=2, line=1.5, cex=0.65)
  } 
  if(tiss=="Kidney_Aged"){
    text="rat kidney aged adult"
  }
  
  text(text, x=rat.start+(rat.end-rat.start)/50, y=2000, cex=0.85, adj=c(0,1))
}

#####################################################################################

## expdiv ortho

expdiv.ortho=expdiv[ortho.neighbours.ortho,"ExpressionDivergence"]
expdiv.spec=expdiv[ortho.neighbours.specific,"ExpressionDivergence"]

d.ortho=density(expdiv.ortho, na.rm=T, bw=0.05)
d.spec=density(expdiv.spec, na.rm=T, bw=0.05)

xlim=range(c(d.ortho$x, d.spec$x))
ylim=range(c(d.ortho$y, d.spec$y))


par(mar=c(3.5, 3.5, 1.5, 0.1)) ## for annotation plots

plot(d.ortho$x, d.ortho$y, col="black", type="l", xlim=xlim, ylim=ylim, axes=F, xlab="", ylab="")
lines(d.spec$x, d.spec$y, col="red", type="l")

pval=wilcox.test(expdiv.ortho, expdiv.spec)$p.value
print(pval)

axis(side=2, cex.axis=0.85, mgp=c(3,0.75,0))
mtext("density", side=2, line=2, cex=0.65)
axis(side=1, cex.axis=0.85, mgp=c(3,0.5,0))
mtext("raw expression divergence", side=1, line=2, cex=0.65)


legend("topright", c("close to ortho lnc", "close to specific lnc"), lty=1, col=c("black", "red"), inset=c(0.01,-0.1), xpd=NA, cex=0.95, bty="n") 

mtext("B", side=3, line=0.5, at=-0.45, font=2, cex=0.9)


#####################################################################################


## expdiv ortho

expdiv.ortho=expdiv[ortho.neighbours.ortho,"ResidualExpressionDivergence"]
expdiv.spec=expdiv[ortho.neighbours.specific,"ResidualExpressionDivergence"]

d.ortho=density(expdiv.ortho, na.rm=T, bw=0.05)
d.spec=density(expdiv.spec, na.rm=T, bw=0.05)

xlim=range(c(d.ortho$x, d.spec$x))
ylim=range(c(d.ortho$y, d.spec$y))


pval=wilcox.test(expdiv.ortho, expdiv.spec)$p.value
print(pval)



par(mar=c(3.5, 3.5, 1.5, 0.1)) ## for annotation plots

plot(d.ortho$x, d.ortho$y, col="black", type="l", xlim=xlim, ylim=ylim, axes=F, xlab="", ylab="")
lines(d.spec$x, d.spec$y, col="red", type="l")

axis(side=2, cex.axis=0.85, mgp=c(3,0.75,0))
mtext("density", side=2, line=2, cex=0.65)
axis(side=1, cex.axis=0.85, mgp=c(3,0.5,0))
mtext("residual expression divergence", side=1, line=2, cex=0.65)


mtext("C", side=3, line=0.5, at=-0.65, font=2, cex=0.9)


#####################################################################################


dev.off()

#####################################################################################
