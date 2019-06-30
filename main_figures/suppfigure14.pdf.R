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

  
  spliced=allinfo.mouse$GeneID[which(allinfo.mouse$NbExons>1)]
  bidirpc=allinfo.mouse$GeneID[which(!is.na(allinfo.mouse$BidirectionalPromoterProteinCoding))]
  bidirother=setdiff(allinfo.mouse$GeneID[which(!is.na(allinfo.mouse$BidirectionalPromoter))], bidirpc)
  unidir=setdiff(allinfo.mouse$GeneID, c(bidirpc, bidirother))
  enhancers=allinfo.mouse$GeneID[which(allinfo.mouse$OverlapEncodeEnhancer=="Yes")]
  antisense=allinfo.mouse$GeneID[which(allinfo.mouse$RegionPCAntisense!="intergenic")]
  intergenic=allinfo.mouse$GeneID[which(allinfo.mouse$RegionPCAntisense=="intergenic")]

  expdiv$ID.Mouse=unlist(lapply(expdiv$ID, function(x) unlist(strsplit(x, split="_"))[1]))
  rownames(expdiv)=expdiv$ID.Mouse

  
  pc=intersect(pc, rownames(expdiv))
  lnc=intersect(lnc, rownames(expdiv))

  lm1=lm(expdiv$ExpressionDivergence~log2(expdiv$MeanTPM+1))
  expdiv$ResidualExpressionDivergence=lm1$residuals
  
  
  prepare=FALSE
}


######################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure14.pdf", sep=""), width=6.85, height=6.25)

#############################################################################

m=matrix(rep(NA, 18*8), nrow=18)


for(i in 1:7){
  m[i,]=c(rep(1,4), rep(2,4))
}


for(i in 8:12){
  m[i,]=2+c(rep(1, 2), rep(3, 2), rep(5, 2), rep(7,2))
}


for(i in 13:17){
  m[i,]=2+c(rep(2, 2), rep(4, 2), rep(6, 2), rep(8,2))
}

m[18,]=2+rep(9,8)


layout(m)

#############################################################################

## raw expression divergence by promoter type


ylim=c(0,1)
xlim=c(0.6,9.9)

xpos=c(1,2, 3.5,4.5, 6,7, 8.5,9.5)

par(mar=c(3.1,3.1,2.2,0.5))

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

boxplot(expdiv[intersect(pc, unidir), "ExpressionDivergence"], boxwex=0.7, col="indianred", at=xpos[1], add=T, axes=F, outline=F, notch=T)
boxplot(expdiv[intersect(lnc, unidir), "ExpressionDivergence"], boxwex=0.7, col="steelblue", at=xpos[2], add=T, axes=F, outline=F, notch=T)

boxplot(expdiv[intersect(pc, bidirpc), "ExpressionDivergence"], boxwex=0.7, col="indianred", at=xpos[3], add=T, axes=F, outline=F, notch=T)
boxplot(expdiv[intersect(lnc, bidirpc), "ExpressionDivergence"], boxwex=0.7, col="steelblue", at=xpos[4], add=T, axes=F, outline=F, notch=T)

boxplot(expdiv[intersect(pc, bidirother), "ExpressionDivergence"], boxwex=0.7, col="indianred", at=xpos[5], add=T, axes=F, outline=F, notch=T)
boxplot(expdiv[intersect(lnc, bidirother), "ExpressionDivergence"], boxwex=0.7, col="steelblue", at=xpos[6], add=T, axes=F, outline=F, notch=T)

boxplot(expdiv[intersect(pc, enhancers), "ExpressionDivergence"], boxwex=0.7, col="indianred", at=xpos[7], add=T, axes=F, outline=F, notch=T)
boxplot(expdiv[intersect(lnc, enhancers), "ExpressionDivergence"], boxwex=0.7, col="steelblue", at=xpos[8], add=T, axes=F, outline=F, notch=T)


axis(side=1, at=xpos, labels=rep("",length(xpos)))
axis(side=2, cex.axis=0.85, mgp=c(3,0.5,0))
mtext("raw expression divergence", side=2, line=2, cex=0.65)

mtext("unidirectional", side=1, at=mean(xpos[1:2]), line=1.25, cex=0.65)
mtext("bidirectional,\nprotein-coding", side=1, at=mean(xpos[3:4]), line=2.25, cex=0.65)

mtext("bidirectional,\nother", side=1, at=mean(xpos[5:6]), line=2.25, cex=0.65)
mtext("enhancer\noverlap", side=1, at=mean(xpos[7:8]), line=2.25, cex=0.65)


mtext("A", side=3, at=-0.8, font=2, cex=0.9, line=0.75)

#############################################################################

## residual expression divergence by promoter type

ylim=c(-0.3,0.8)
xlim=c(0.6,9.9)

xpos=c(1,2, 3.5,4.5, 6,7, 8.5,9.5)

par(mar=c(3.1,3.1,2.2,0.5))

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)


boxplot(expdiv[intersect(pc, unidir), "ResidualExpressionDivergence"], boxwex=0.7, col="indianred", at=xpos[1], add=T, axes=F, outline=F, notch=T)
boxplot(expdiv[intersect(lnc, unidir), "ResidualExpressionDivergence"], boxwex=0.7, col="steelblue", at=xpos[2], add=T, axes=F, outline=F, notch=T)

boxplot(expdiv[intersect(pc, bidirpc), "ResidualExpressionDivergence"], boxwex=0.7, col="indianred", at=xpos[3], add=T, axes=F, outline=F, notch=T)
boxplot(expdiv[intersect(lnc, bidirpc), "ResidualExpressionDivergence"], boxwex=0.7, col="steelblue", at=xpos[4], add=T, axes=F, outline=F, notch=T)

boxplot(expdiv[intersect(pc, bidirother), "ResidualExpressionDivergence"], boxwex=0.7, col="indianred", at=xpos[5], add=T, axes=F, outline=F, notch=T)
boxplot(expdiv[intersect(lnc, bidirother), "ResidualExpressionDivergence"], boxwex=0.7, col="steelblue", at=xpos[6], add=T, axes=F, outline=F, notch=T)

boxplot(expdiv[intersect(pc, enhancers), "ResidualExpressionDivergence"], boxwex=0.7, col="indianred", at=xpos[7], add=T, axes=F, outline=F, notch=T)
boxplot(expdiv[intersect(lnc, enhancers), "ResidualExpressionDivergence"], boxwex=0.7, col="steelblue", at=xpos[8], add=T, axes=F, outline=F, notch=T)


axis(side=1, at=xpos, labels=rep("",length(xpos)))
axis(side=2, cex.axis=0.85, mgp=c(3,0.5,0))
mtext("residual expression divergence", side=2, line=2, cex=0.65)

mtext("unidirectional", side=1, at=mean(xpos[1:2]), line=1.25, cex=0.65)
mtext("bidirectional,\nprotein-coding", side=1, at=mean(xpos[3:4]), line=2.25, cex=0.65)

mtext("bidirectional,\nother", side=1, at=mean(xpos[5:6]), line=2.25, cex=0.65)
mtext("enhancer\noverlap", side=1, at=mean(xpos[7:8]), line=2.25, cex=0.65)


legend("topleft", legend=c("protein-coding", "lncRNA"), fill=c("indianred", "steelblue"), bty="n", inset=c(-0.01,-0.01), cex=0.95, horiz=T, xpd=NA)

mtext("B", side=3, at=-0.8, font=2, cex=0.9, line=0.75)

#############################################################################

## examples of genes with highly divergent expression

samples=kronecker(tissue.order, age.order, paste, sep="_")
xpos=1:20
names(xpos)=samples
xpos=setdiff(xpos, 16)
names(xpos)=setdiff(samples, "Testis_EarlyEmbryo")


sptissage=colnames(avgexp.mr)
sp=unlist(lapply(sptissage, function(x) unlist(strsplit(x, split="_"))[1]))
tissage=unlist(lapply(sptissage, function(x) paste(unlist(strsplit(x, split="_"))[2:3], collapse="_")))
tiss=unlist(lapply(sptissage, function(x) unlist(strsplit(x, split="_"))[2]))

residuals.pc=expdiv$ResidualExpressionDivergence[which(expdiv$GeneType=="protein_coding")]
names(residuals.pc)=expdiv$ID[which(expdiv$GeneType=="protein_coding")]

residuals.lnc=expdiv$ResidualExpressionDivergence[which(expdiv$GeneType=="lncRNA")]
names(residuals.lnc)=expdiv$ID[which(expdiv$GeneType=="lncRNA")]

residuals.pc=sort(residuals.pc, decreasing=T)
residuals.lnc=sort(residuals.lnc, decreasing=T)

## 2 pc examples, 2 lnc examples

for(i in 1:4){
  if(i<3){
    example=names(residuals.pc)[i]
    mouse=unlist(strsplit(example, split="_"))[1]
    name=allinfo.mouse[mouse, "GeneName"]
    locus=paste(allinfo.mouse[mouse, "Chr"],":",allinfo.mouse[mouse, "Start"],"-",allinfo.mouse[mouse, "End"],":",allinfo.mouse[mouse, "Strand"], sep="")
  } else{
    example=names(residuals.lnc)[i-2]
    mouse=unlist(strsplit(example, split="_"))[1]
    name=allinfo.mouse[mouse, "GeneName"]
    locus=paste(allinfo.mouse[mouse, "Chr"],":",allinfo.mouse[mouse, "Start"],"-",allinfo.mouse[mouse, "End"],":",allinfo.mouse[mouse, "Strand"], sep="")
  }
  
  ## actual plot
  
  smallx=0.1
  
  for(this.sp in c("Mouse", "Rat")){
    ylim=c(0, max(log2(as.numeric(avgexp.mr[example,which(sp==this.sp)])+1)))
    xlim=c(0.5,20.5)
    

    if(i==1){
      if(this.sp=="Mouse"){
        par(mar=c(0.1, 3.5, 2.1, 0.0))
      } else{
        par(mar=c(2.1, 3.5, 0.1, 0.0))
      }
    } else{
      if(this.sp=="Mouse"){
        par(mar=c(0.1, 2, 2.1, 1.5))
      } else{
        par(mar=c(2.1, 2, 0.1, 1.5))
      }
    }
    
    plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, xlab="", ylab="")
    
    rect(xpos[tissage[which(sp==this.sp)]]-smallx, 0, xpos[tissage[which(sp==this.sp)]]+smallx, log2(as.numeric(avgexp.mr[example,which(sp==this.sp)])+1), col=col.tissues[tiss[which(sp==this.sp)]], border=col.tissues[tiss[which(sp==this.sp)]])
    if(this.sp=="Rat"){
      axis(side=1, at=1:20, labels=rep("",20), mgp=c(3,0.5,0), cex.axis=0.7)
      mtext(rep(1:5, 4), side=1, at=1:20, cex=0.5, line=0.5)
    }
    
    axis(side=2, cex.axis=0.85, mgp=c(3, 0.5, 0))

   
    if(this.sp=="Mouse"){
      if(!is.na(name)){
        mtext(name, side=3, font=3, cex=0.65)
      } else{
        mtext(mouse, side=3, font=3, cex=0.65)
      }
    }
   
    if(i==1){
      text(tolower(this.sp), x=1, adj=c(0,0), y=ylim[2]+diff(ylim)/100, cex=1, xpd=NA)

      if(this.sp=="Mouse"){
        mtext("avg. expression level (log2-transformed TPM)", side=2, line=2, at=0-ylim[2]/15, cex=0.65)

        mtext("C", side=3, at=-6.25, font=2, cex=0.8, line=0.8)
      } 
       
    }

    if(this.sp=="Rat"){
        mtext("developmental stage", side=1, line=1.5, cex=0.65)
      }
  }
}


######################################################################

par(mar=c(0,0,0,0))
plot(1, type="n", xlab="", ylab="", axes=F)
  
legend("topleft", legend=shortname.tiss, fill=col.tissues, border=col.tissues, inset=c(0.05,0.035), cex=0.95, horiz=T, bty="n") 

######################################################################

dev.off()

######################################################################
