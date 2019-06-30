#####################################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE
  
  minTPM=1
  maxFDR=0.01
}


###################################################################################

if(load==TRUE){
  for(sp in c("Mouse", "Rat", "Chicken")){
    load(paste("RData/data.annotations.", sp, ".RData", sep=""))
    assign(paste("pc", tolower(sp), sep="."), pc)
    assign(paste("lnc", tolower(sp), sep="."), lnc)
    assign(paste("allinfo", tolower(sp), sep="."), allinfo)
    
    load(paste("RData/data.expression.",sp,".RData", sep=""))
    assign(paste("normtpm", tolower(sp), sep="."), normtpm)

    load(paste("RData/data.tpm.stats.",sp,".RData", sep=""))
    assign(paste("stats", tolower(sp), sep="."), stats) 

    if(sp!="Chicken"){
      load(paste("RData/data.diffexp.",sp,".RData", sep=""))
      assign(paste("de.allreads.", tolower(sp), sep=""), de.global.allreads)
      assign(paste("de.resampled.", tolower(sp), sep=""), de.global.resampled)
    }
    
  }

  load=FALSE
}

###################################################################################

if(prepare==TRUE){
  
  table.maxstage=list()
  
  foldchange=list("pc"=list("Mouse"=list(), "Rat"=list()),"lnc"=list("Mouse"=list(), "Rat"=list()))
  

  for(sp in c("Mouse", "Rat", "Chicken")){
    avgtpm=get(paste("stats", tolower(sp), sep="."))
    
    this.normtpm=get(paste("normtpm", tolower(sp), sep="."))

    this.info=get(paste("allinfo", tolower(sp), sep="."))
    this.highexp=this.info$GeneID[which(this.info$MaxTPM>=minTPM)]
    
    this.pc=get(paste("pc", tolower(sp), sep="."))
    this.lnc=get(paste("lnc", tolower(sp), sep="."))

    ## take only genes above noise expression level
    this.pc=intersect(this.pc, this.highexp)
    this.lnc=intersect(this.lnc, this.highexp)
    
    maxexp=log2(apply(this.normtpm, 1, function(x) max(x, na.rm=T))+1)

    assign(paste("maxexp", tolower(sp), sep="."), maxexp)
    assign(paste("pc", tolower(sp), sep="."), this.pc)
    assign(paste("lnc", tolower(sp), sep="."), this.lnc)

    if(sp%in%c("Mouse", "Rat")){
      this.de=get(paste("de.allreads", tolower(sp), sep="."))

      
      table.maxstage[[sp]]=list()

      for(tiss in tissue.order){
        this.age.order=age.order

        if(tiss=="Testis"){
          this.age.order=setdiff(age.order, "EarlyEmbryo")
        }
        
        this.tpm=avgtpm[,paste("MeanTPM.",tiss, "_", this.age.order,sep="")]
        colnames(this.tpm)=this.age.order
        this.maxstage=apply(this.tpm,1, function(x) this.age.order[which.max(x)])
        names(this.maxstage)=rownames(this.tpm)

        this.maxtpm=apply(this.tpm,1,max)
        this.mintpm=apply(this.tpm,1,min)

        this.ratio=100*(this.maxtpm-this.mintpm)/this.maxtpm
        names(this.ratio)=rownames(this.tpm)
        
        this.signif=this.de$GeneID[which(this.de[,paste("FDR", tiss, sep=".")]<maxFDR)]
        this.signif.pc=intersect(this.pc, this.signif)
        this.signif.lnc=intersect(this.lnc, this.signif)


        foldchange[["pc"]][[sp]][[tiss]]=this.ratio[this.signif.pc]
        foldchange[["lnc"]][[sp]][[tiss]]=this.ratio[this.signif.lnc]

        this.table.pc=as.numeric(table(factor(this.maxstage[this.signif.pc], levels=age.order)))
        names(this.table.pc)=age.order

        this.table.pc=100*this.table.pc/sum(this.table.pc)

        this.table.lnc=as.numeric(table(factor(this.maxstage[this.signif.lnc], levels=age.order)))
        names(this.table.lnc)=age.order

        this.table.lnc=100*this.table.lnc/sum(this.table.lnc)

        table.maxstage[[sp]][[tiss]]=list("pc"=this.table.pc, "lnc"=this.table.lnc)
      }
    }
    
  }

  prepare=FALSE
}

###################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure9.pdf", sep=""), width=6.85, height=2.75)

###################################################################################

## layout

m=matrix(rep(NA,10*11), nrow=10)

for(i in 1:10){
  m[i,]=c(rep(1,5), rep(2,6))
}

layout(m)

###################################################################################

par(mar=c(3,2.75,2.25,1))

plot(1, type="n", axes=F, xlab="", ylab="", ylim=c(0, 105), xlim=c(0.5,4.5), xaxs="i")

xpos=c(1,2,3,4)
names(xpos)=tissue.order
pos.sp=c()

for(tiss in tissue.order){
  this.pos=xpos[tiss]+c(-0.3,-0.15, 0.15, 0.3)
  pos.sp=c(pos.sp, xpos[tiss]+c(-0.225,0.225))
  boxplot(foldchange[["pc"]][["Mouse"]][[tiss]], foldchange[["lnc"]][["Mouse"]][[tiss]], foldchange[["pc"]][["Rat"]][[tiss]], foldchange[["lnc"]][["Rat"]][[tiss]], add=T, at=this.pos, xlab="", ylab="", axes=F, outline=F, boxwex=0.12, col=c(col.tissues[tiss],"white"), border=c("black",col.tissues[tiss]), notch=T)
}

axis(side=2, mgp=c(3,0.5,0), cex.axis=0.85, at=pretty(c(0,100)))
mtext("relative expression change (% max)", side=2, line=1.5, cex=0.65)

axis(side=1, mgp=c(3,0.5,0), cex.axis=0.75, at=pos.sp, labels=rep("", length(pos.sp)))
mtext(rep(c("mouse","rat"), length(pos.sp)/2), line=0.5, side=1, cex=0.55, at=pos.sp)


legend("topleft", legend=shortname.tiss, fill=col.tissues, border=col.tissues, bty="o", box.col="white", inset=c(0.01,-0.125), cex=0.95, bg="white", horiz=T, xpd=NA)


mtext("A", side=3, line=1, font=2, cex=0.95, at=0.1)


###################################################################################


## maxstage for diffexp pc and lnc, for rat only

this.table.maxstage=table.maxstage[["Rat"]]

addspace=0.5
xpos=c(1:5, 6:10+addspace, 11:15+2*addspace, 16:20+3*addspace)
names(xpos)=kronecker(tissue.order, age.order, paste, sep="_")

ylim=c(0,60)
xlim=c(min(xpos)-0.5,max(xpos)+0.5)

width=0.2
space=0.04

par(mar=c(3.15, 2.75, 2.5, 0.25))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", xaxs="i")

for(tiss in tissue.order){
  for(age in age.order){
    this.x=xpos[paste(tiss,age, sep="_")]
    this.prop.pc=this.table.maxstage[[tiss]][["pc"]][age]
    this.prop.lnc=this.table.maxstage[[tiss]][["lnc"]][age]

    if(this.prop.pc>0 & this.prop.lnc>0){
      rect(this.x-width-space, 0, this.x-space, this.prop.pc, col=col.tissues[tiss], border=col.tissues[tiss])
      rect(this.x+space, 0, this.x+width+space, this.prop.lnc, col="white", border=col.tissues[tiss])
    }
  }
}

axis(side=2,  at=seq(from=0, to=60, by=20), cex.axis=0.9, mgp=c(3, 0.5, 0))
axis(side=1, at=xpos, labels=rep("",20), mgp=c(3,0.5,0))
mtext(rep(1:5, 4), at=xpos, line=0.5, cex=0.6, side=1)
mtext("% differentially expressed genes", side=2, line=2, cex=0.65)
mtext("developmental stage with maximum expression, rat genes", side=1, line=1.8, cex=0.65)

abline(v=c(mean(xpos[5:6]), mean(xpos[10:11]), mean(xpos[15:16])), lty=2, col="gray40")

legend("topright", legend=c("protein-coding", "lncRNAs"), fill=c("black", "white"), bty="n", cex=0.95, inset=c(-0.01,-0.125), xpd=NA)


mtext("B", side=3, line=0.5, at=-1.5, font=2, cex=0.95)

###################################################################################


dev.off()

###################################################################################
