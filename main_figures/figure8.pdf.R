##############################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE

  library(ape)

  maxFDR=0.01
}

##########################################################################

if(load==TRUE){

  print("loading data")

  load("RData/data.expression.ortho.RData")
  
  tpm.ortholnc=exportho.mr[which(exportho.mr$GeneType=="lncRNA"), which(!colnames(exportho.mr)%in%c("ID", "GeneType"))]

  idmouse.ortholnc=unlist(lapply(rownames(tpm.ortholnc), function(x) unlist(strsplit(x, split="_"))[1]))
  idrat.ortholnc=unlist(lapply(rownames(tpm.ortholnc), function(x) unlist(strsplit(x, split="_"))[2]))
  id.ortholnc=rownames(tpm.ortholnc)

  
  load("RData/data.diffexp.Mouse.RData")
  de.mouse=de.global.allreads

  load("RData/data.diffexp.Rat.RData")
  de.rat=de.global.allreads

  load(paste("RData/data.tpm.stats.Mouse.RData", sep=""))
  avgtpm.mouse=stats

  load(paste("RData/data.tpm.stats.Rat.RData", sep=""))
  avgtpm.rat=stats

  ## kmeans for differentially expressed lncRNAs, both species
  
  kmeans=list()
  
  for(tiss in c("Brain","Kidney", "Liver", "Testis")){
    this.kmeans=read.table(paste(pathDatasets, "SupplementaryDataset4/KmeansClusters_MouseRatOrtho_DiffExp_",tiss,"_LncRNAs.txt", sep=""), h=T,stringsAsFactors=F)

    kmeans[[tiss]]=this.kmeans
  }
  
  load=FALSE
}

##############################################################################

if(prepare==TRUE){

  shared.diffexp.lnc=list()

  for(tiss in tissue.order){
    this.age.order=age.order
    
    if(tiss=="Testis"){
      this.age.order=setdiff(age.order, "EarlyEmbryo")
    }

    this.tpm.mouse=avgtpm.mouse[,paste("MeanTPM.",tiss, "_", this.age.order,sep="")]
    colnames(this.tpm.mouse)=this.age.order
    this.maxstage.mouse=apply(this.tpm.mouse,1, function(x) this.age.order[which.max(x)])
    names(this.maxstage.mouse)=rownames(this.tpm.mouse)

    this.tpm.rat=avgtpm.rat[,paste("MeanTPM.",tiss, "_", this.age.order,sep="")]
    colnames(this.tpm.rat)=this.age.order
    this.maxstage.rat=apply(this.tpm.rat,1, function(x) this.age.order[which.max(x)])
    names(this.maxstage.rat)=rownames(this.tpm.rat)
        
    signif.mouse=de.mouse$GeneID[which(de.mouse[,paste("FDR", tiss, sep=".")]<maxFDR)]
    signif.rat=de.rat$GeneID[which(de.rat[,paste("FDR", tiss, sep=".")]<maxFDR)]
    
    signif.both=which(idmouse.ortholnc%in%signif.mouse & idrat.ortholnc%in%signif.rat)

    this.id.mouse=idmouse.ortholnc[signif.both]
    this.id.rat=idrat.ortholnc[signif.both]

    allres=data.frame("ID.Mouse"=this.id.mouse,"ID.Rat"=this.id.rat, "MaxStage.Mouse"=this.maxstage.mouse[this.id.mouse], "MaxStage.Rat"=this.maxstage.rat[this.id.rat], stringsAsFactors=F)

    shared.diffexp.lnc[[tiss]]=allres
    
  }
  
  
  prepare=FALSE
 }

###################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "Figure8.pdf", sep=""), width=6.85, height=10.5)

###################################################################################

m=matrix(rep(NA,16*10), nrow=16)

for(i in 1:4){
  m[i,]=c(rep(1,10))
}

for(i in 5:7){
  m[i,]=c(rep(2:6, each=2))
}

for(i in 8:10){
  m[i,]=c(rep(7:11, each=2))
}

for(i in 11:13){
  m[i,]=c(rep(12:16, each=2))
}

for(i in 14:16){
  m[i,]=c(rep(17:21, each=2))
}

layout(m)

par(oma=c(1,2.1,0,0))

###################################################################################

addspace=0.5
xpos=c(1:5, 6:10+addspace, 11:15+2*addspace, 16:20+3*addspace)
names(xpos)=kronecker(tissue.order, age.order, paste, sep="_")


ylim=c(0,100)
xlim=c(min(xpos)-0.7,max(xpos)+2.75)

width=0.3
space=0.04

par(mar=c(4.1,3.1,2.1,3.1))
     
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", xaxs="i")

for(tiss in tissue.order){
  this.de=shared.diffexp.lnc[[tiss]]

  for(age in age.order){
    this.de.age=as.numeric(table(factor(this.de$MaxStage.Rat[which(this.de$MaxStage.Mouse==age)], levels=age.order)))
    names(this.de.age)=levels(age.order)

    this.de.age=100*this.de.age/sum(this.de.age)

    this.xpos=xpos[paste(tiss,age,sep="_")]
    
    rect(this.xpos-width, c(0, cumsum(this.de.age)[1:4]), this.xpos+width, cumsum(this.de.age), col=col.ages, border=col.ages)

    n=length(which(this.de$MaxStage.Mouse==age))
    
    mtext(n, at=this.xpos, side=1, line=2.5, cex=0.67)
  }
  
}

mtext("N =", side=1, line=2.5, at=-1, cex=0.67)

axis(side=1, at=xpos, labels=rep("",length(xpos)))
axis(side=2, cex.axis=0.85, mgp=c(3,0.5,0))

mtext(rep(1:5, 4), at=xpos, side=1, line=0.5, cex=0.65)
mtext("stage with maximum expression in mouse", side=1, line=1.5, cex=0.65, at=mean(xpos))
mtext("% max. expression in stage, rat", side=2, cex=0.7, line=2.5, at=50, outer=F, xpd=NA)

mtext("A", side=3, line=1, cex=0.9, at=-1.5, font=2)

mtext("brain", side=3, at=mean(xpos[1:5]), cex=0.75, line=0)
mtext("kidney", side=3, at=mean(xpos[6:10]), cex=0.75, line=0)
mtext("liver", side=3, at=mean(xpos[11:15]), cex=0.75, line=0)
mtext("testes", side=3, at=mean(xpos[17:20]), cex=0.75, line=0)


legend("topright", legend=c(1:5), fill=col.ages, border=col.ages, bty="n", inset=c(0.01,-0.01), cex=0.95, horiz=F, xpd=NA)


###################################################################################

labels=c("B","C", "D", "E")
names(labels)=c("Brain", "Kidney", "Liver", "Testis")

mean.xpos=apply(matrix(c(1,1.25, 2,2.25, 3, 3.25, 4, 4.25, 5, 5.25),nrow=2, byrow=F),2,mean)

for(tiss in c("Brain", "Kidney", "Liver", "Testis")){

  if(tiss=="Testis"){
    xpos=c(2,2.25, 3, 3.25, 4, 4.25, 5, 5.25)
  } else{
    xpos=c(1,1.25, 2,2.25, 3, 3.25, 4, 4.25, 5, 5.25)
  }
  
  this.kmeans=kmeans[[tiss]]
  ## order genes by cluster
  
  for(i in 1:length(unique(this.kmeans$ClusterIndex))){
    
    clust=i
    
    this.genes=which(this.kmeans$ClusterIndex==clust)
    relexp=this.kmeans[this.genes,-c(1:3)]
    
    meanexp=apply(relexp,2,mean)
    minexp=apply(relexp,2,min)
    maxexp=apply(relexp,2,max)

    
    ## par(mar=c(2.1,2.5,1.75, 0.75))

    if(i==1){
      par(mar=c(1.1,0.5,1.75, 0.5))
    } else{
      par(mar=c(1.1,0.5,1.75, 0.5))
    }
    
    plot(1, xlim=c(0.9, 5.35), ylim=c(0,1), axes=F, type="n")
    
    mtext(1:5, at=mean.xpos, side=1, line=0.4, cex=0.6)

    odds=seq(from=1, to=length(xpos), by=2)
    evens=seq(from=2, to=length(xpos), by=2)
        
    for(g in this.genes){
      lines(xpos[odds], as.numeric(relexp[g,odds]), col="gray80", lwd=0.5)
      lines(xpos[evens], as.numeric(relexp[g,evens]), col="gray80", lwd=0.5)
    }
        
    points(xpos[odds], meanexp[odds], pch= pch.allsp["mouse"], bg=col.tissues[tiss], type="b")
    points(xpos[evens], meanexp[evens], pch= pch.allsp["rat"], bg=col.tissues[tiss], type="b")
    
    box()
        
    axis(side=1, at=mean.xpos, labels=rep("",length(mean.xpos)))

    if(i==1){
      axis(side=2, cex.axis=0.85, mgp=c(3, 0.5, 0), at=c(0,0.5,1))
      mtext(labels[tiss], side=3, at=-0.2, line=0.5, font=2, cex=0.9)
    }

    if(i==1){
      mtext("relative expression", side=2, line=1.5, cex=0.65, outer=F, xpd=NA)
    }
    
    if(i==3 ){     
      mtext("developmental stage", side=1, cex=0.65, line=1.1)
    }
    

    smallx=diff(xpos)[1]
    
    if(which.max(meanexp)%in%c(1:3)){
      text(paste("N=",length(this.genes),sep=""), x=xpos[length(xpos)-1]-smallx, y=0.95, cex=0.95)
    } else{
      text(paste("N=",length(this.genes),sep=""), x=xpos[length(xpos)-1]-smallx, y=0.05, cex=0.95)
    }
    
  }
}

###################################################################################

## legend for tissue and species
plot(1, type="n",xlab="", ylab="", axes=F)


legend("topleft", legend=c("mouse","rat"), pch=pch.allsp[c("mouse","rat")], bty="n", cex=0.95, inset=c(0.05, -0.05),  xpd=NA)


legend("topleft", legend=shortname.tiss, fill=col.tissues, border=col.tissues, bty="n", inset=c(0.05,0.35), cex=0.95,  xpd=NA)

###################################################################################



###################################################################################

dev.off()

###################################################################################

