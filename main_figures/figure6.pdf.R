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
  load("RData/data.pca.RData")
  load("RData/data.expression.conservation.spearman.RData")
  load("RData/data.expression.clustering.RData")

  load=FALSE
}

##############################################################################

if(prepare==TRUE){
  samples=rownames(pca.pc$co)
  species=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))
  tissues=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[2]))
  ages=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[3]))
  ages=unlist(lapply(ages, function(x) substr(x, 1, nchar(x)-1)))
  
 
  explained.lnc=round(100*pca.lnc$eig/sum(pca.lnc$eig), digits=2)

  prepare=FALSE
}

##############################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "Figure6.pdf", sep=""), width=6.85, height=3.5)

#############################################################################

## layout

m=matrix(rep(NA, 11*12), nrow=11)

for(i in 1:10){
  m[i,]=c(rep(1,5), rep(2,2), rep(3,4), rep(4,1))
}


for(i in 11){
  m[i,]=c(rep(5,12))
}


layout(m)

#############################################################################


## first plot: PCA for lncRNAs

par(mar=c(3.75, 3.5, 2.5, 0.5))

this.pch=pch.allsp[tolower(species)]
this.col=col.tissage[paste(tissues, ages, sep="_")]

plot(pca.lnc$co[,1], pca.lnc$co[,2], pch=this.pch, bg=this.col, col="black", cex=1.35, xlab="", ylab="", axes=F)

axis(side=1, cex.axis=0.85, mgp=c(3,0.5, 0))
axis(side=2, cex.axis=0.85, mgp=c(3,0.65, 0))
box()
mtext(paste("PC1 (",explained.lnc[1],"% explained variance)", sep=""), side=1, line=1.65, cex=0.65)
mtext(paste("PC2 (",explained.lnc[2],"% explained variance)", sep=""), side=2, line=2.1, cex=0.65)
## mtext("PCA, protein-coding genes", side=3, cex=0.7, font=2, line=0.25)

legend("topright", legend=c("mouse","rat"), pch=pch.allsp[1:2], bty="n", cex=0.95)

mtext("A", side=3, at=-0.975, line=1.1, font=2, cex=0.9)

#############################################################################


## dendrogram species lnc
par(tck=NA)

par(mar=c(2.89,0.5,1.2,0))

plot(tree.lnc.tpm, direction="rightwards",show.tip.label=F,edge.width=0.75)

mtext("B", side=3, at=0, line=-0.05, font=2, cex=0.9)

###############################################################################


## heatmap for lncRNAs, Spearman, hierarchical clustering
zlim=range(c(as.numeric(cormat.lnc.tpm)))

par(mar=c(3.5, 0.1, 1.8, 1.5))
image(cormat.lnc.tpm[sample.order.lnc.tpm,sample.order.lnc.tpm],col=terrain.colors(75), zlim=zlim,axes=F, xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")

nbsamples=length(sample.order.lnc.tpm)
width=1/nbsamples
xstart=seq(from=0, to=1-width, length=nbsamples)
this.col=col.tissage[paste(tissues, ages, sep="_")[sample.order.lnc.tpm]]

rect(xstart, -0.06, xstart+width, -0.02,  col=this.col, border=NA, xpd=NA)

this.col=col.allsp[tolower(species[sample.order.lnc.tpm])]

rect(1.02, xstart, 1.06, xstart+width, col=this.col, border=NA, xpd=NA)

mtext("species", side=4, line=1.15, cex=0.65)

mtext("organ / developmental stage", side=1, line=1.25, cex=0.65)


###############################################################################

## legend

par(mar=c(0,0.5,2.5,0))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")

## this.tissage=kronecker(c("br", "kd", "lv", "ts"), 1:5, paste, sep=" ")
## this.tissage=setdiff(this.tissage, "ts 1")

## xpos=0.1
## xwidth=0.15
## ywidth=0.01

## ypos=seq(from=1, to=0.4, length=20)
## ypos=ypos[-16]

## rect(xpos, ypos, xpos+xwidth, ypos+ywidth, col=col.tissage, border=NA, xpd=NA)
## text(unlist(lapply(this.tissage, function(x) unlist(strsplit(x, split=" "))[2])), x=xpos+xwidth*1.8, y=ypos, adj=c(0.5,0), cex=0.8, xpd=NA)
## text("br", x=xpos+4*xwidth, y=mean(ypos[1:5]), cex=0.9, adj=c(0.5,0))
## text("kd", x=xpos+4*xwidth, y=mean(ypos[6:10]), cex=0.9, adj=c(0.5,0))
## text("lv", x=xpos+4*xwidth, y=mean(ypos[11:15]), cex=0.9, adj=c(0.5,0))
## text("ts", x=xpos+4*xwidth, y=mean(ypos[16:19]), cex=0.9, adj=c(0.5,0))


xpos=-0.1
xwidth=0.2
ywidth=0.02
ypos=seq(from=0.25, to=0.18, length=2)
rect(xpos, ypos, xpos+xwidth, ypos+ywidth, col=col.allsp, border=NA, xpd=NA)
text(c("mouse", "rat"), x=xpos+2*xwidth, y=ypos,  cex=0.95, adj=c(0,0))

###############################################################################


## legend for correlations

par(mar=c(1.1,36.25,0.35,9))

nbcol=75
breaks=nbcol+1

min.val=zlim[1]
max.val=zlim[2]

z <- seq(min.val, max.val, length = nbcol)
xax=pretty(z)
xax=xax[which(xax>=min(z) & xax<=max(z))]
image(x=z, z = matrix(z, ncol = 1), col = terrain.colors(nbcol), zlim=zlim, xlim=range(xax), xaxt="n" ,yaxt="n")

par(tck=-0.35)
axis(side=1, at = xax, labels = xax,cex.axis=0.75,mgp=c(3,0.1,0))

mtext("Spearman's\ncorrelation",side=4, at=0.0, cex=0.65,line=0.5, las=2)

##############################################################################

this.tissage=kronecker(c("br", "kd", "lv", "ts"), 1:5, paste, sep=" ")
this.tissage=setdiff(this.tissage, "ts 1")

ypos=0.5
ywidth=0.75
xwidth=0.1

xpos=seq(from=-5, to=-0.5, length=20)
xpos=xpos[-16]

rect(xpos,  ypos,  xpos+xwidth,ypos+ywidth, col=col.tissage, border=NA, xpd=NA)
text(unlist(lapply(this.tissage, function(x) unlist(strsplit(x, split=" "))[2])), x=xpos+xwidth/2, y=ypos-ywidth/3, adj=c(0.5,1), cex=0.85, xpd=NA)
text("br", y=ypos-2.5*ywidth, x=mean(xpos[1:5]), cex=0.95, adj=c(0.5,1), xpd=NA)
text("kd", y=ypos-2.5*ywidth, x=mean(xpos[6:10]), cex=0.95, adj=c(0.5,1), xpd=NA)
text("lv", y=ypos-2.5*ywidth, x=mean(xpos[11:15]), cex=0.95, adj=c(0.5,1), xpd=NA)
text("ts", y=ypos-2.5*ywidth, x=mean(xpos[16:19]), cex=0.95, adj=c(0.5,1), xpd=NA)

dev.off()

##############################################################################
