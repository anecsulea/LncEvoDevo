###################################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE
  
  pch.allsp=c(21, 23, 24)
  names(pch.allsp)=c("mouse", "rat", "chicken")

  col.allsp=c("black", "gray40", "gray80")
  names(col.allsp)=c("mouse", "rat", "chicken")
  

  library(ape)
}

###################################################################################

if(load==TRUE){

  load("RData/data.pca.mrc.RData")
  load("RData/data.expression.clustering.mrc.RData")
  load("RData/data.sample.quality.RData")

  load=FALSE
}

###################################################################################

if(prepare==TRUE){

  samples=rownames(pca.pc$co)
  species=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))
  tissues=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[2]))
  ages=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[3]))
  ages=unlist(lapply(ages, function(x) substr(x, 1, nchar(x)-1)))
  
  explained.pc=round(100*pca.pc$eig/sum(pca.pc$eig), digits=2)
  explained.devTF=round(100*pca.pc.devTF$eig/sum(pca.pc.devTF$eig), digits=2)

  print(all(rownames(pca.pc.devTF)==rownames(pca.pc)))

  deg=deg[rownames(pca.pc$co),]
   
  prepare=FALSE
}

###################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure6.pdf", sep=""), width=6.85, height=9.5)

###################################################################################

## layout

m=matrix(rep(NA,23*24), nrow=23)

for(i in 1:8){
  m[i,]=c(rep(1,11), rep(2,3), rep(3,10))
}

for(i in 9){
  m[i,]=c(rep(4,16), rep(5,8))
}

for(i in 10:16){
  m[i,]=c(rep(6,6), rep(7,6), rep(8,6), rep(9,6))
}


for(i in 17:23){
  m[i,]=c(rep(10,6), rep(11, 6),rep(12,6), rep(13, 6))
}

layout(m)

###################################################################################


## first plot: PCA for protein-coding genes

par(mar=c(3.75, 3.5, 2.5, 0.5))

this.pch=pch.allsp[tolower(species)]
this.col=col.tissage[paste(tissues, ages, sep="_")]

plot(pca.pc.devTF$co[,1], pca.pc.devTF$co[,2], pch=this.pch, bg=this.col, col="black", cex=1.35, xlab="", ylab="", axes=F)

axis(side=1, cex.axis=0.95, mgp=c(3,0.5, 0))
axis(side=2, cex.axis=0.95, mgp=c(3,0.65, 0))
box()
mtext(paste("PC1 (",explained.devTF[1],"% explained variance)", sep=""), side=1, line=2, cex=0.7)
mtext(paste("PC2 (",explained.devTF[2],"% explained variance)", sep=""), side=2, line=2.25, cex=0.7)

legend("topright", legend=c("mouse","rat", "chicken"), pch=pch.allsp, bty="n", cex=1)

mtext("A", side=3, at=-0.95, line=1.25, font=2, cex=0.95)

#############################################################################

## dendrogram species pc

par(tck=NA)

par(mar=c(1.75,1,2.7,0))

plot(tree.devTF.tpm, direction="rightwards",show.tip.label=F,edge.width=0.75)

mtext("B", side=3, at=0.01, line=1, font=2, cex=0.95)

###############################################################################

## heatmap for protein-coding genes, Spearman, hierarchical clustering

par(mar=c(2.5, 0.1, 3.5, 2.25))

zlim=range(c(as.numeric(cormat.devTF.tpm)))

image(cormat.devTF.tpm[sample.order.devTF.tpm,sample.order.devTF.tpm],col=terrain.colors(75), zlim=zlim,axes=F, xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")

nbsamples=length(sample.order.devTF.tpm)
width=1/nbsamples
xstart=seq(from=0, to=1-width, length=nbsamples)
this.col=col.tissage[paste(tissues, ages, sep="_")[sample.order.devTF.tpm]]
rect(xstart, -0.06, xstart+width, -0.02,  col=this.col, border=NA, xpd=NA)

## species 

this.col=col.allsp[tolower(species[sample.order.devTF.tpm])]

rect(1.02, xstart, 1.07, xstart+width, col=this.col, border=NA, xpd=NA)

mtext("species", side=4, line=1.25, cex=0.7)

mtext("organ / developmental stage", side=1, line=1.5, cex=0.7)

legend("topright", legend=c("mouse", "rat", "chicken"), fill=col.allsp, border=col.allsp, bty="n", inset=c(0,-0.1), horiz=T, xpd=NA)

###############################################################################

## legend

par(mar=c(0,0.5,1.5,5.5))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")

ypos=0.7
xpos=seq(from=0.05, to=1, length=20)
xpos=xpos[-16]

ywidth=0.2
xwidth=0.025


this.tissage=kronecker(c("br", "kd", "lv", "ts"), 1:5, paste, sep=" ")
this.tissage=setdiff(this.tissage, "ts 1")

rect(xpos, ypos,  xpos+xwidth, ypos+ywidth, col=col.tissage, border=NA, xpd=NA)
text(unlist(lapply(this.tissage, function(x) unlist(strsplit(x, split=" "))[2])), y=ypos+ywidth*2.75, x=xpos+xwidth/2, adj=c(0.5,0.5), cex=0.9, xpd=NA)
text("brain", y=ypos+4.25*ywidth, x=mean(xpos[1:5]), cex=0.95, adj=c(0.5,0), xpd=NA)
text("kidney", y=ypos+4.25*ywidth, x=mean(xpos[6:10]), cex=0.95, adj=c(0.5,0), xpd=NA)
text("liver", y=ypos+4.25*ywidth, x=mean(xpos[11:15]), cex=0.95, adj=c(0.5,0), xpd=NA)
text("testes", y=ypos+4.25*ywidth, x=mean(xpos[16:19]), cex=0.95, adj=c(0.5,0), xpd=NA)

###############################################################################

## legend for correlations

par(mar=c(1.25,3,1,6.5))

nbcol=75
breaks=nbcol+1
## zlim was previously defined
min.val=zlim[1]
max.val=zlim[2]

z <- seq(min.val, max.val, length = nbcol)
xax=pretty(z)
xax=xax[which(xax>=min(z) & xax<=max(z))]
image(x=z, z = matrix(z, ncol = 1), col = terrain.colors(nbcol), zlim=zlim, xlim=range(xax), xaxt="n" ,yaxt="n")

par(tck=-0.5)
axis(side=1, at = xax, labels = xax, cex.axis=0.8, mgp=c(3,0.5,0))

mtext("Spearman's\ncorrelation",side=4, at=0.0, cex=0.6,line=0.5, las=2)

par(tck=-0.05)

###############################################################################

## correlation between developmental stages and PC1/PC2

pc1=pca.pc.devTF$co[,1]
pc2=pca.pc.devTF$co[,2]

for(tiss in tissue.order){

  if(tiss=="Brain" | tiss=="Testis"){
    this.pc=pc1
    i=1
  } else{
    this.pc=pc2
    i=2
  }
  
  ylim=range(this.pc[which(tissues==tiss)])
  xlim=c(0.75, 5.25)
  
  xpos=1:5
  names(xpos)=age.order
  
  par(mar=c(4.1, 3.1, 4.1, 0.5))
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)
  
  for(age in age.order){
    this.sp=species[which(ages==age & tissues==tiss)]
    this.x=rep(xpos[age], length(which(ages==age & tissues==tiss)))
    
    points(this.x, this.pc[which(ages==age & tissues==tiss)], pch=pch.allsp[tolower(this.sp)], col="black", bg=col.tissues[tiss], cex=1.35)
}

  box()

  axis(side=1, at=1:5, cex.axis=0.85, labels=rep("", 5), mgp=c(3,0.5,0))
  axis(side=2, cex.axis=0.85, mgp=c(3,0.75,0))
  
  mtext(paste("coordinates on PC", i, sep=""), side=2, line=2.1, cex=0.7)
  
  mtext(tolower(tiss), side=3, line=0.15, cex=0.75)
  
  mtext(1:5, side=1, line=0.5, cex=0.7, at=1:5)

  if(tiss=="Brain"){
    mtext("C", side=3, line=1.8, font=2, cex=0.95, at=-0.65)
  }
}

mtext("coordinates on PCA axes", side=3, oma=T, at=-7, cex=0.75, line=2)
mtext("developmental stages", side=1, oma=T, at=-7, cex=0.75, line=2)

###############################################################################

## correlation between developmental stages and PC1/PC2, after correcting for coverage bias

pc1=pca.pc.devTF$co[,1]
pc2=pca.pc.devTF$co[,2]

for(tiss in tissue.order){

  if(tiss=="Brain" | tiss=="Testis"){
    this.pc=lm(pc1~deg$CoverageBias)$residuals
    i=1
  } else{
    this.pc=lm(pc2~deg$CoverageBias)$residuals
    i=2
  }  
  
  ylim=range(this.pc[which(tissues==tiss)])
  xlim=c(0.75, 5.25)
  
  xpos=1:5
  names(xpos)=age.order
  
  par(mar=c(4.1, 3.1, 4.1, 0.5))
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)
  
  for(age in age.order){
    this.sp=species[which(ages==age & tissues==tiss)]
    this.x=rep(xpos[age], length(which(ages==age & tissues==tiss)))
    
    points(this.x, this.pc[which(ages==age & tissues==tiss)], pch=pch.allsp[tolower(this.sp)], col="black", bg=col.tissues[tiss], cex=1.35)
}

  box()

  axis(side=1, at=1:5, cex.axis=0.85, labels=rep("", 5), mgp=c(3,0.75,0))
  axis(side=2, cex.axis=0.85, mgp=c(3,0.75,0))
  
  mtext(paste("coordinates on PC", i, ", residuals", sep=""), side=2, line=2.1, cex=0.7)
  
  mtext(tolower(tiss), side=3, line=0.15, cex=0.75)

  mtext(1:5, side=1, line=0.5, cex=0.7, at=1:5)

  if(tiss=="Brain"){
    mtext("D", side=3, line=1.8, font=2, cex=0.95, at=-0.65)
  }
}


mtext("coordinates on PCA axes, after correcting for 3' read coverage bias", side=3, oma=T, at=-7, cex=0.75, line=2)
mtext("developmental stages", side=1, oma=T, at=-7, cex=0.75, line=2)

###############################################################################

dev.off()

###############################################################################

