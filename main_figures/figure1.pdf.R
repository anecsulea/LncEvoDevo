###################################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE
}

library("imager") ## to include graphics
library(ape)

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

  load("RData/data.expression.ortho.RData")
  avgexp.mr=avgexp.mr[,setdiff(colnames(avgexp.mr), c("GeneID", "GeneType"))]

  
  load("RData/data.celltype.markers.RData")

  ## PCA and clustering
  
  load("RData/data.pca.mrc.RData")
  load("RData/data.expression.clustering.mrc.RData")

  
  load=F
}

###################################################################################

if(prepare==TRUE){
  
  mouse.img=load.image("images/mouse.png")
  rat.img=load.image("images/rat.jpg")
  chicken.img=load.image("images/chicken.jpg")
  allsp.img=load.image("images/Silhouettes.png")
  
  ## data for PCA

  samples=rownames(pca.pc$co)
  species=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))
  tissues=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[2]))
  ages=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[3]))
  ages=unlist(lapply(ages, function(x) substr(x, 1, nchar(x)-1)))
  
  explained.pc=round(100*pca.pc$eig/sum(pca.pc$eig), digits=2)
  explained.devTF=round(100*pca.pc.devTF$eig/sum(pca.pc.devTF$eig), digits=2)
  
  prepare=F
}

###################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "Figure1.pdf", sep=""), width=6.85, height=5.75)

###################################################################################

## layout

m=matrix(rep(NA,10*24), nrow=10)

for(i in 1:3){
  m[i,]=c(rep(1,2), rep(2,22))
}

for(i in 4:9){
  m[i,]=c(rep(3,11), rep(4,3), rep(5,10))
}

for(i in 10){
  m[i,]=c(rep(6,16), rep(7,8))
}

layout(m)


#####################################################################################

par(mar=c(0., 0.2, 1.15, 1.15))
plot(allsp.img, axes=F)

mtext("A", side=3, at=155, font=2, cex=0.95, line=0)

###################################################################################

## plot: developmental stages

par(mar=c(1.75, 0.0, 0.75, 0.1))

plot(1, type="n", xlab="", ylab="", xlim=c(-0.5,6), ylim=c(0.2,3), axes=F, xaxs="i", yaxs="i")

ypos=2.75
ym=2.15
yr=1.4
yc=0.7

x1=1.35
x2=2.6
x3=3.8
x4=4.5
x5=5.5
xmax=5.75
smally=0.1
largey=0.2

segments(0, ypos, x1, ypos, lwd=1)
segments(x1, ypos, x2, ypos, lty=2)
segments(x2, ypos, x3, ypos)
segments(x3, ypos, x4, ypos, lty=2)
segments(x4, ypos, x5, ypos)

xaxs=c(0.675, 0.875, 1.15, 3.25, 5) ## developmental stages


segments(c(0,1), ypos-largey, c(0,1), ypos+largey, col="black") ## conception and birth
segments(xaxs, ypos-smally, xaxs, ypos+smally, col="blue", lwd=1.25)

arrows(x5, ypos, xmax, ypos, length=0.1, angle=50, lwd=1)

xaxslarge=c(-0.5, 0.6, 1.65) ## where to place the labels

segments(xaxslarge[1:3]+0.12, ym+0.2, xaxs[1:3], ypos-smally-0.02, col="gray", lty=2)
## text for developmental stages

text(x=xaxslarge+c(0.05,0,0), y=ym, labels=c("E13.5", "E17.5", "P1-2"), adj=0, cex=0.95)
text(x=xaxslarge+c(0.055,0,0), y=yr, labels=c("E15", "E18.5", "P1-2"), adj=0, cex=0.95)
text(x=xaxslarge[1:2]+c(0.05,0), y=yc, labels=c("HH31", "HH36"), adj=0, cex=0.95)

text(x=xaxs[4:5], y=ym, labels=c("8-10 weeks", "12-24 months"), cex=0.95)
text(x=xaxs[4:5], y=yr, labels=c("8-10 weeks", "12-24 months"), cex=0.95)


## text(x=0, y=3.1, "conception", col="black", xpd=NA, cex=0.95)
## text(x=1, y=3.1, "birth", col="black", xpd=NA, cex=0.95)
## text(x=1, y=3.12, "0", col="black", xpd=NA, cex=0.95)

text(as.character(1:5), x=xaxs, y=3.05, col="black", xpd=NA, cex=0.95)

## tissue colors

xt=0.25
yt=0.27

## tissues for first stage
text(c("br", "kd", "lv"), col=col.tissues[c("Brain", "Kidney", "Liver")], x=c(xaxslarge[1], xaxslarge[1]+xt, xaxslarge[1]+2*xt), y=ym-yt, adj=0.5, cex=1.1, xpd=NA)
text(c("br", "kd", "lv"), col=col.tissues[c("Brain", "Kidney", "Liver")], x=c(xaxslarge[1], xaxslarge[1]+xt, xaxslarge[1]+2*xt), y=yr-yt, adj=0.5, cex=1.1, xpd=NA)
text(c("br", "kd", "lv"), col=col.tissues[c("Brain", "Kidney", "Liver")], x=c(xaxslarge[1], xaxslarge[1]+xt, xaxslarge[1]+2*xt), y=yc-yt, adj=0.5, cex=1.1, xpd=NA)

## second stage

xt2=0.22
xt3=0.1
text(c("br", "kd", "lv", "ts"), col=col.tissues[c("Brain", "Kidney", "Liver", "Testis")], x=xt3+c(xaxslarge[2]-xt2, xaxslarge[2], xaxslarge[2]+xt2, xaxslarge[2]+2*xt2), y=ym-yt, adj=0.5, cex=1.1)
text(c("br", "kd", "lv", "ts"), col=col.tissues[c("Brain", "Kidney", "Liver", "Testis")], x=xt3+c(xaxslarge[2]-xt2, xaxslarge[2], xaxslarge[2]+xt2, xaxslarge[2]+2*xt2), y=yr-yt, adj=0.5, cex=1.1)
text(c("br", "kd", "lv"), col=col.tissues[c("Brain", "Kidney", "Liver")], x=xt3+c(xaxslarge[2]-xt2, xaxslarge[2], xaxslarge[2]+xt2), y=yc-yt, adj=0.5, cex=1.1)



## third stage

xt2=0.22
xt3=0.1
text(c("br", "kd", "lv", "ts"), col=col.tissues[c("Brain", "Kidney", "Liver", "Testis")], x=xt3+c(xaxslarge[3]-xt2, xaxslarge[3], xaxslarge[3]+xt2, xaxslarge[3]+2*xt2), y=ym-yt, adj=0.5, cex=1.1)
text(c("br", "kd", "lv", "ts"), col=col.tissues[c("Brain", "Kidney", "Liver", "Testis")], x=xt3+c(xaxslarge[3]-xt2, xaxslarge[3], xaxslarge[3]+xt2, xaxslarge[3]+2*xt2), y=yr-yt, adj=0.5, cex=1.1)

## fourth and fifth stage

for(i in 4:5){
  xt2=0.12
  xt3=0.05
  text(c("br", "kd", "lv", "ts"), col=col.tissues[c("Brain", "Kidney", "Liver", "Testis")], x=xt3+c(xaxs[i]-3*xt2, xaxs[i]-xt2, xaxs[i]+xt2, xaxs[i]+3*xt2), y=ym-yt, adj=0.5, cex=1.1)
  text(c("br", "kd", "lv", "ts"), col=col.tissues[c("Brain", "Kidney", "Liver", "Testis")], x=xt3+c(xaxs[i]-3*xt2, xaxs[i]-xt2, xaxs[i]+xt2, xaxs[i]+3*xt2), y=yr-yt, adj=0.5, cex=1.1)
}
 
## labels for stages
segments(xaxslarge[1]-0.5*xt, yc-2*yt, xaxslarge[1]+2.5*xt, yc-2*yt, xpd=NA, lwd=1, col="black")
text("mid-stage embryo", x=xaxslarge[1]+xt,  y=yc-2.65*yt, cex=0.95, adj=0.5, xpd=NA)
## text("mid-stage", x=xaxslarge[1]+xt,  y=yc-2.6*yt, cex=0.9, adj=0.5, xpd=NA)
## text("embryo", x=xaxslarge[1]+xt,  y=yc-3.6*yt, cex=0.9, adj=0.5, xpd=NA)

segments(xaxslarge[2]-0.75*xt, yc-2*yt, xaxslarge[2]+2.5*xt, yc-2*yt, xpd=NA, lwd=1, col="black")
text("late embryo", x=xaxslarge[2]+xt*0.9,  y=yc-2.65*yt, cex=0.95, adj=0.5, xpd=NA)
## text("late", x=xaxslarge[2]+xt,  y=yc-2.6*yt, cex=0.9, adj=0.5, xpd=NA)
## text("embryo", x=xaxslarge[2]+xt,  y=yc-3.6*yt, cex=0.9, adj=0.5, xpd=NA)

segments(xaxslarge[3]-0.75*xt, yc-2*yt, xaxslarge[3]+2.6*xt, yc-2*yt, xpd=NA, lwd=1, col="black")
text("newborn", x=xaxslarge[3]+xt*0.8,  y=yc-2.65*yt, cex=0.95, adj=0.5, xpd=NA)

segments(xaxs[4]-3.5*xt2, yc-2*yt, xaxs[4]+4*xt2, yc-2*yt, xpd=NA, lwd=1, col="black")
text("young adult", x=xaxs[4],  y=yc-2.65*yt, cex=0.95, adj=0.5, xpd=NA)

segments(xaxs[5]-3.5*xt2, yc-2*yt, xaxs[5]+4*xt2, yc-2*yt, xpd=NA, lwd=1, col="black")
text("aged adult", x=xaxs[5],  y=yc-2.65*yt, cex=0.95, adj=0.5, xpd=NA)

###################################################################################

## first plot: PCA for protein-coding genes

par(mar=c(2.75, 3.5, 4.5, 0.5))

this.pch=pch.allsp[tolower(species)]
this.col=col.tissage[paste(tissues, ages, sep="_")]

plot(pca.pc$co[,1], pca.pc$co[,2], pch=this.pch, bg=this.col, col="black", cex=1.35, xlab="", ylab="", axes=F)

axis(side=1, cex.axis=0.95, mgp=c(3,0.5, 0))
axis(side=2, cex.axis=0.95, mgp=c(3,0.65, 0))
box()
mtext(paste("PC1 (",explained.pc[1],"% explained variance)", sep=""), side=1, line=2, cex=0.7)
mtext(paste("PC2 (",explained.pc[2],"% explained variance)", sep=""), side=2, line=2.25, cex=0.7)


legend("topright", legend=c("mouse","rat", "chicken"), pch=pch.allsp, bty="n", cex=1)

mtext("B", side=3, at=-0.971, line=1, font=2, cex=0.95)

#############################################################################

## dendrogram species pc

par(tck=NA)

par(mar=c(1.7,1,3.7,0))

plot(tree.pc.tpm, direction="rightwards",show.tip.label=F,edge.width=0.75)

mtext("C", side=3, at=0.01, line=0.1, font=2, cex=0.95)

###############################################################################


## heatmap for protein-coding genes, Spearman, hierarchical clustering

par(mar=c(2.5, 0.1, 4.5, 2.25))
zlim=range(c(as.numeric(cormat.pc.tpm)))

image(cormat.pc.tpm[sample.order.pc.tpm,sample.order.pc.tpm],col=terrain.colors(75), zlim=zlim,axes=F, xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")

nbsamples=length(sample.order.pc.tpm)
width=1/nbsamples
xstart=seq(from=0, to=1-width, length=nbsamples)

this.col=col.tissage[paste(tissues, ages, sep="_")[sample.order.pc.tpm]]

rect(xstart, -0.06, xstart+width, -0.02,  col=this.col, border=NA, xpd=NA)

## species 

this.col=col.allsp[tolower(species[sample.order.pc.tpm])]

rect(1.02, xstart, 1.07, xstart+width, col=this.col, border=NA, xpd=NA)

mtext("species", side=4, line=1.25, cex=0.7)

mtext("organ / developmental stage", side=1, line=1.5, cex=0.7)


legend("topright", legend=c("mouse", "rat", "chicken"), fill=col.allsp, border=col.allsp, bty="n", inset=c(0,-0.1), horiz=T, xpd=NA)

###############################################################################

## legend

par(mar=c(0,0.5,2.5,5.5))
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

par(mar=c(2.25,2.5,1,6.5))

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
axis(side=1, at = xax, labels = xax, cex.axis=0.85, mgp=c(3,0.5,0))

mtext("Spearman's\ncorrelation",side=4, at=0.0, cex=0.6,line=0.5, las=2)

###############################################################################

dev.off()

#######################################################################################

