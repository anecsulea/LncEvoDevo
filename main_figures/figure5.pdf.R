##########################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE

}

##########################################################################

if(load==TRUE){
  ## only from mouse perspective

  sp="Mouse"
  load(paste("RData/data.annotations.", sp, ".RData", sep=""))

  load(paste("RData/data.tpm.stats.",sp,".RData", sep=""))
  avgtpm=stats

  load(paste("RData/data.expression.",sp,".RData",sep=""))

  load("RData/data.projection.stats.Mouse.Rat.RData")
  proj.mouse.rat=projstats

  load("RData/data.projection.stats.Mouse.Chicken.RData")
  proj.mouse.chicken=projstats

  ## ortho families
  load("RData/data.ortho.families.RData")

  load("RData/data.expression.ortho.RData")

  load=FALSE
}  

##########################################################################

if(prepare==TRUE){


  ## max age and tissue for mouse, based on averaged TPM values across replicates
  
  maxsample=avgtpm[,"MaxSample"]
  maxtissue=unlist(lapply(maxsample, function(x) unlist(strsplit(x, split="_"))[1]))
  maxage=unlist(lapply(maxsample, function(x) unlist(strsplit(x, split="_"))[2]))
  
  names(maxsample)=rownames(avgtpm)
  names(maxtissue)=rownames(avgtpm)
  names(maxage)=rownames(avgtpm)

  ## ortho
  
  minratio=0
  minreads=10

  ## vs rat
  
  projected.unfiltered.rat=proj.mouse.rat$GeneID[which((proj.mouse.rat$Length.UnfilteredProjections/proj.mouse.rat$TotalExonicLength)>minratio)]

  projected.step1.rat=proj.mouse.rat$GeneID[which((proj.mouse.rat$Length.FilteredProjectionsStep1/proj.mouse.rat$TotalExonicLength)>minratio)]

  projected.step2.rat=proj.mouse.rat$GeneID[which((proj.mouse.rat$Length.FilteredProjectionsStep2/proj.mouse.rat$TotalExonicLength)>minratio)]

  projected.transcribed.rat=proj.mouse.rat$GeneID[which((proj.mouse.rat$Length.FilteredProjectionsStep2/proj.mouse.rat$TotalExonicLength)>minratio & proj.mouse.rat$TotalReadCount.TargetSpecies.AllSamples>minreads)]

  ortho.rat=orthofam$ID.Mouse[which(!is.na(orthofam$ID.Rat))]

  ortholnc.rat=orthofam$ID.Mouse[which(!is.na(orthofam$ID.Rat) & orthofam$GeneType.Rat=="lncRNA")]
  orthopc.rat=orthofam$ID.Mouse[which(!is.na(orthofam$ID.Rat) & orthofam$GeneType.Rat=="protein_coding")]

  ## vs chicken
  
  projected.unfiltered.chicken=proj.mouse.chicken$GeneID[which((proj.mouse.chicken$Length.UnfilteredProjections/proj.mouse.chicken$TotalExonicLength)>minratio)]
  
  projected.step1.chicken=proj.mouse.chicken$GeneID[which((proj.mouse.chicken$Length.FilteredProjectionsStep1/proj.mouse.chicken$TotalExonicLength)>minratio)]

  projected.step2.chicken=proj.mouse.chicken$GeneID[which((proj.mouse.chicken$Length.FilteredProjectionsStep2/proj.mouse.chicken$TotalExonicLength)>minratio)]

  projected.transcribed.chicken=proj.mouse.chicken$GeneID[which((proj.mouse.chicken$Length.FilteredProjectionsStep2/proj.mouse.chicken$TotalExonicLength)>minratio & proj.mouse.chicken$TotalReadCount.TargetSpecies.AllSamples>minreads)]

  ortho.chicken=orthofam$ID.Mouse[which(!is.na(orthofam$ID.Chicken))]

  ortholnc.chicken=orthofam$ID.Mouse[which(!is.na(orthofam$ID.Chicken) & orthofam$GeneType.Chicken=="lncRNA")]
  orthopc.chicken=orthofam$ID.Mouse[which(!is.na(orthofam$ID.Chicken) & orthofam$GeneType.Chicken=="protein_coding")]

  ## genes with no ortho
  
  lnc.noortho=setdiff(lnc, c(ortholnc.rat, ortholnc.chicken))
  pc.noortho=setdiff(pc, c(orthopc.rat, orthopc.chicken))


  ortholnc.all=orthofam[which((!is.na(orthofam$ID.Chicken) & orthofam$GeneType.Chicken=="lncRNA" & (!is.na(orthofam$ID.Rat)) & orthofam$GeneType.Rat=="lncRNA" & (!is.na(orthofam$ID.Mouse)) & orthofam$GeneType.Mouse=="lncRNA")),]
  

  prepare=TRUE

}

##########################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "Figure5.pdf", sep=""), width=6.85, height=3.5)

m=matrix(rep(NA, 5*20), nrow=5)

for(i in 1:5){
  m[i,]=c(rep(1,8), rep(2, 6), rep(3, 6))
}

layout(m)

###################################################################################

## barplot: stats for mouse vs. rat

nb.lnc=c(length(lnc), length(intersect(lnc, projected.step2.rat)), length(intersect(lnc, projected.transcribed.rat)), length(intersect(lnc, ortho.rat)), length(intersect(lnc, ortholnc.rat)), length(intersect(lnc, projected.step2.chicken)), length(intersect(lnc, projected.transcribed.chicken)), length(intersect(lnc, ortho.chicken)), length(intersect(lnc, ortholnc.chicken)))

nb.pc=c(length(pc), length(intersect(pc, projected.step2.rat)), length(intersect(pc, projected.transcribed.rat)), length(intersect(pc, ortho.rat)), length(intersect(pc, orthopc.rat)), length(intersect(pc, projected.step2.chicken)), length(intersect(pc, projected.transcribed.chicken)), length(intersect(pc, ortho.chicken)), length(intersect(pc, orthopc.chicken)))

width=0.12
xpos=c(1, 2:5+0.5, 6:9+1)
xlnc=0.32
xlim=c(0.65, 10.75)
ylim=c(0, 22000)


par(mar=c(8, 3.1, 2.5, 0.1))
plot(1, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=F, xaxs="i", yaxs="i")

rect(xpos-width, 0, xpos+width, nb.pc, col="indianred", border="indianred")
rect(xpos-width+xlnc, 0, xpos+width+xlnc, nb.lnc, col="steelblue", border="steelblue")

yax=pretty(ylim/1000)
yax=yax[which(yax>=(ylim[1]/1000) & yax<=(ylim[2]/1000))]
axis(side=2, cex.axis=0.85, mgp=c(3, 0.75, 0), at=yax*1000, labels=paste(yax,"K",sep=""))

mtext("number of genes", side=2, cex=0.65, line=2)
axis(side=1, at=xpos+xlnc/2, labels=rep("", 9))


segments(xpos[2], -7400, xpos[5], -7400, xpd=NA)
mtext("mouse vs. rat", side=1, line=6, at=mean(xpos[2:5]), cex=0.65)


segments(xpos[6], -7400, xpos[9], -7400, xpd=NA)
mtext("mouse vs. chicken", side=1, line=6, at=mean(xpos[6:9]), cex=0.65)

smally=150
text(labels=nb.lnc, x=xpos+xlnc+width/3, y=nb.lnc+smally, srt=90, cex=0.85, adj=c(0,0.5))

mtext(c("all genes", rep(c("conserved\nsequence", "transcribed", "putative ortho", "confirmed\northo"), 2)), side=1, las=2, line=0.65, cex=0.55, at=xpos+xlnc/2)

legend("topright", c("protein-coding", "lncRNA"), fill=c("indianred","steelblue"),  border=c("indianred","steelblue"), inset=c(0.05, -0.05), bty="n", cex=0.95, xpd=NA)

mtext("A", side=3, at=-0.75, line=1, font=2, cex=0.9)

##########################################################################

## table for tissue 
tab.pc.noortho=as.numeric(table(factor(maxtissue[pc], levels=tissue.order)))
tab.pc.ortho=as.numeric(table(factor(maxtissue[intersect(pc, ortho.rat)], levels=tissue.order)))
tab.pc.ortho.chicken=as.numeric(table(factor(maxtissue[intersect(pc, ortho.chicken)], levels=tissue.order)))

tab.lnc.noortho=as.numeric(table(factor(maxtissue[lnc], levels=tissue.order)))
tab.lnc.ortho=as.numeric(table(factor(maxtissue[intersect(lnc, ortho.rat)], levels=tissue.order)))
tab.lnc.ortho.chicken=as.numeric(table(factor(maxtissue[intersect(lnc, ortho.chicken)], levels=tissue.order)))


m=matrix(c(tab.pc.noortho, tab.lnc.noortho,  tab.pc.ortho, tab.lnc.ortho,  tab.pc.ortho.chicken, tab.lnc.ortho.chicken), nrow=6, byrow=6)
m=m/apply(m,1,sum)

xpos=c(1, 1.35, 2, 2.35, 3, 3.35)

width=0.15

par(mar=c(4.1,3,2.5,0.1))

plot(1, xlim=c(0.8,4.25), ylim=c(0,100), axes=F, xlab="", ylab="", type="n")

for(i in 1:6){
  ypos=100*cumsum(m[i,])
  rect(xpos[i]-width, c(0, ypos[-length(ypos)]), xpos[i]+width, ypos, col=col.tissues, border=NA)
}

axis(side=2, cex.axis=0.85, mgp=c(3,0.75,0))
axis(side=1, at=xpos, labels=rep("",6))

mtext(rep(c("pc", "lnc"), 3), side=1, line=0.5, cex=0.65, at=xpos)
mtext(rep("ortho",2), side=1, line=1.6, at=c(mean(xpos[3:4]),mean(xpos[5:6])), cex=0.65)
mtext(c("no ortho","rat", "chicken"), side=1, line=c(2, 2.4, 2.4), cex=0.65, at=c(mean(xpos[1:2]),mean(xpos[3:4]),mean(xpos[5:6])))
mtext("% genes w. max. exp. in organ", side=2, cex=0.65, line=1.9)

legend("topright", legend=c("ts", "lv", "kd", "br"), border=col.tissues[4:1], fill=col.tissues[4:1], horiz=F, bty="n", xpd=NA, inset=c(-0.01, 0.01), cex=0.95)

mtext("B", side=3, at=-0.05, line=0.75, font=2, cex=0.9)


##########################################################################

## stage with maximum expression

tab.pc.noortho=as.numeric(table(factor(maxage[pc.noortho], levels=age.order)))
tab.pc.ortho=as.numeric(table(factor(maxage[intersect(pc, ortho.rat)], levels=age.order)))
tab.pc.ortho.chicken=as.numeric(table(factor(maxage[intersect(pc, ortho.chicken)], levels=age.order)))

tab.lnc.noortho=as.numeric(table(factor(maxage[lnc.noortho], levels=age.order)))
tab.lnc.ortho=as.numeric(table(factor(maxage[intersect(lnc, ortho.rat)], levels=age.order)))
tab.lnc.ortho.chicken=as.numeric(table(factor(maxage[intersect(lnc, ortho.chicken)], levels=age.order)))


m=matrix(c(tab.pc.noortho, tab.lnc.noortho,  tab.pc.ortho, tab.lnc.ortho,  tab.pc.ortho.chicken, tab.lnc.ortho.chicken), nrow=6, byrow=6)
m=m/apply(m,1,sum)

xpos=c(1, 1.35, 2, 2.35, 3, 3.35)

width=0.15

par(mar=c(4.1,3,2.5,0.1))

plot(1, xlim=c(0.8,4.25), ylim=c(0,100), axes=F, xlab="", ylab="", type="n")

for(i in 1:6){
  ypos=100*cumsum(m[i,])
  rect(xpos[i]-width, c(0, ypos[-length(ypos)]), xpos[i]+width, ypos, col=col.ages, border=NA)
}

axis(side=2, cex.axis=0.85, mgp=c(3,0.75,0))
axis(side=1, at=xpos, labels=rep("",6))

mtext("% genes w. max. exp. in stage", side=2, cex=0.65, line=1.9)

mtext(rep(c("pc", "lnc"), 3), side=1, line=0.5, cex=0.65, at=xpos)

mtext(rep("ortho",2), side=1, line=1.6, at=c(mean(xpos[3:4]),mean(xpos[5:6])), cex=0.65)
mtext(c("no ortho","rat", "chicken"), side=1, line=c(2, 2.4, 2.4), cex=0.65, at=c(mean(xpos[1:2]),mean(xpos[3:4]),mean(xpos[5:6])))

legend("topright", legend=5:1, border=col.ages[5:1], fill=col.ages[5:1], horiz=F, bty="n", xpd=NA, cex=0.95, inset=c(-0.01,0.01))

mtext("C", side=3, at=-0.05, line=0.75, font=2, cex=0.9)


##########################################################################


dev.off()

##########################################################################
