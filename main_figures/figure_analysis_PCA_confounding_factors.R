###############################################################################

source("parameters.R")

pathHisat=paste(path, "results/hisat/", sep="")
pathDegradation=paste(path, "results/RNA_degradation/", sep="")

###############################################################################

source("parameters.R")


###############################################################################

deg=read.table(paste(pathDegradation, "Stats_RIN_3bias.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

deg$Sample.ID=paste(deg$Species, deg$Sample.ID, sep="_")

rownames(deg)=deg$Sample.ID

colnames(deg)=c("Species", "SampleID", "NbReads", "RIN", "CoverageBias")

###############################################################################

load("RData/data.pca.mrc.RData")

###############################################################################

samples=rownames(pca.pc$co)

deg=deg[samples,]

###############################################################################

species=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))
tissues=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[2]))
ages=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[3]))
ages=unlist(lapply(ages, function(x) substr(x, 1, nchar(x)-1)))

this.pch=pch.allsp[tolower(species)]
this.col=col.tissage[paste(tissues, ages, sep="_")]

###############################################################################

pdf(file=paste(pathFigures, "PCA_ConfoundingFactors.pdf", sep=""), width=5, height=6)

m=matrix(1:6, nrow=3, byrow=T)
layout(m)

###############################################################################

## PCA vs. RIN

par(mar=c(3.1, 4.1, 2.1, 2.1))

plot(pca.pc$co[,1], deg$RIN, pch=this.pch, bg=this.col, col="black", xlab="", ylab="", axes=F, ylim=c(6.75,10.25))

axis(side=1, cex.axis=0.85, mgp=c(3,0.5, 0))
axis(side=2, cex.axis=0.85, mgp=c(3,0.65, 0))
box()
mtext(paste("PC1", sep=""), side=1, line=2, cex=0.7)
mtext(paste("RIN", sep=""), side=2, line=2.15, cex=0.7)

r2=round(cor(pca.pc$co[,1], deg$RIN)^2, digits=2)
mtext(paste("R2 =", r2), side=3, cex=0.7)


par(mar=c(3.1, 4.1, 2.1, 2.1))

plot(pca.pc$co[,2], deg$RIN, pch=this.pch, bg=this.col, col="black", xlab="", ylab="", axes=F, ylim=c(6.75,10.25))

axis(side=1, cex.axis=0.85, mgp=c(3,0.5, 0))
axis(side=2, cex.axis=0.85, mgp=c(3,0.65, 0))
box()
mtext(paste("PC2", sep=""), side=1, line=2, cex=0.7)
mtext(paste("RIN", sep=""), side=2, line=2.15, cex=0.7)


r2=round(cor(pca.pc$co[,2], deg$RIN)^2, digits=2)
mtext(paste("R2 =", r2), side=3, cex=0.7)

###############################################################################

## PCA vs. 3' bias

par(mar=c(3.1, 4.1, 2.1, 2.1))

plot(pca.pc$co[,1], deg$CoverageBias, pch=this.pch, bg=this.col, col="black", xlab="", ylab="", axes=F)

axis(side=1, cex.axis=0.85, mgp=c(3,0.5, 0))
axis(side=2, cex.axis=0.85, mgp=c(3,0.65, 0))
box()
mtext(paste("PC1", sep=""), side=1, line=2, cex=0.7)
mtext(paste("3' coverage bias", sep=""), side=2, line=2.15, cex=0.7)

r2=round(cor(pca.pc$co[,1], deg$CoverageBias)^2, digits=2)
mtext(paste("R2 =", r2), side=3, cex=0.7)

## PCA vs. 3' bias

par(mar=c(3.1, 4.1, 2.1, 2.1))

plot(pca.pc$co[,2], deg$CoverageBias, pch=this.pch, bg=this.col, col="black", xlab="", ylab="", axes=F)

axis(side=1, cex.axis=0.85, mgp=c(3,0.5, 0))
axis(side=2, cex.axis=0.85, mgp=c(3,0.65, 0))
box()
mtext(paste("PC2", sep=""), side=1, line=2, cex=0.7)
mtext(paste("3' coverage bias", sep=""), side=2, line=2.15, cex=0.7)

r2=round(cor(pca.pc$co[,2], deg$CoverageBias)^2, digits=2)
mtext(paste("R2 =", r2), side=3, cex=0.7)

###############################################################################


## PCA vs. nb unique reads

par(mar=c(3.1, 4.1, 2.1, 2.1))

plot(pca.pc$co[,1], deg$NbReads/1e6, pch=this.pch, bg=this.col, col="black", xlab="", ylab="", axes=F)

axis(side=1, cex.axis=0.85, mgp=c(3,0.5, 0))
axis(side=2, cex.axis=0.85, mgp=c(3,0.65, 0))
box()
mtext(paste("PC1", sep=""), side=1, line=2, cex=0.7)
mtext(paste("nb. million reads", sep=""), side=2, line=2.15, cex=0.7)

r2=round(cor(pca.pc$co[,1], deg$NbReads/1e6)^2, digits=2)
mtext(paste("R2 =", r2), side=3, cex=0.7)



## PCA vs. nb unique reads

par(mar=c(3.1, 4.1, 2.1, 2.1))

plot(pca.pc$co[,2], deg$NbReads/1e6, pch=this.pch, bg=this.col, col="black", xlab="", ylab="", axes=F)

axis(side=1, cex.axis=0.85, mgp=c(3,0.5, 0))
axis(side=2, cex.axis=0.85, mgp=c(3,0.65, 0))
box()
mtext(paste("PC2", sep=""), side=1, line=2, cex=0.7)
mtext(paste("nb. million reads", sep=""), side=2, line=2.15, cex=0.7)

r2=round(cor(pca.pc$co[,2], deg$NbReads/1e6)^2, digits=2)
mtext(paste("R2 =", r2), side=3, cex=0.7)


###############################################################################

dev.off()
    
###############################################################################
