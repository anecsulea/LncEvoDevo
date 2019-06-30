###################################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE
}

###################################################################################

if(load==TRUE){
  for(sp in c("Mouse", "Rat", "Chicken")){
    load(paste("RData/data.annotations.", sp, ".RData", sep=""))
    assign(paste("pc", tolower(sp), sep="."), pc)
    assign(paste("lnc", tolower(sp), sep="."), lnc)
    assign(paste("allinfo", tolower(sp), sep="."), allinfo)
    
    load(paste("RData/data.diffexp.",sp,".RData", sep=""))
    assign(paste("de.consecutive.stages.", tolower(sp), sep=""), de.consecutive.stages)

    load(paste("RData/data.tpm.stats.",sp,".RData", sep=""))
    assign(paste("avgtpm.", tolower(sp), sep=""), stats)
  }

  load=FALSE
}

###################################################################################

if(prepare==TRUE){

  nbuppc=list()
  nbdownpc=list()
  nbuplnc=list()
  nbdownlnc=list()

  maxFDR=0.01
  
  for(sp in c("Mouse", "Rat", "Chicken")){

    this.tissue.order=tissue.order

    if(sp=="Chicken"){
      this.tissue.order=c("Brain","Kidney","Liver")
    }
    
    nbuppc[[sp]]=list()
    nbdownpc[[sp]]=list()
    nbuplnc[[sp]]=list()
    nbdownlnc[[sp]]=list()
    
    this.pc=get(paste("pc", tolower(sp), sep="."))
    this.lnc=get(paste("lnc", tolower(sp), sep="."))
    this.de=get(paste("de.consecutive.stages.", tolower(sp), sep=""))
    this.exp=get(paste("avgtpm.", tolower(sp), sep=""))

    for(tiss in this.tissue.order){
      if(tiss=="Testis"){
        this.ages=c("LateEmbryo", "Newborn", "Adult", "Aged")
      } else{
        if(sp=="Chicken"){
          this.ages=c("EarlyEmbryo", "LateEmbryo")
        } else{
          this.ages=c("EarlyEmbryo", "LateEmbryo", "Newborn", "Adult", "Aged")
        }
      }
      
      nbstages=length(this.ages)

      for(i in 1:(nbstages-1)){
        age1=this.ages[i]
        age2=this.ages[i+1]

        this.signif=this.de$GeneID[which(this.de[,paste("FDR.",tiss,"_",age1,"_",age2,sep="")]<maxFDR)]
        this.signif=intersect(this.signif, rownames(this.exp))
        this.up=this.signif[which(this.exp[this.signif,paste("MeanTPM.",tiss,"_",age2,sep="")]>this.exp[this.signif,paste("MeanTPM.",tiss,"_",age1,sep="")])]
        this.down=this.signif[which(this.exp[this.signif,paste("MeanTPM.",tiss,"_",age1,sep="")]>this.exp[this.signif,paste("MeanTPM.",tiss,"_",age2,sep="")])]

        nbuppc[[sp]][paste(tiss,age1,age2,sep="_")]=length(intersect(this.up, this.pc))
        nbdownpc[[sp]][paste(tiss,age1,age2,sep="_")]=length(intersect(this.down, this.pc))
        
        nbuplnc[[sp]][paste(tiss,age1,age2,sep="_")]=length(intersect(this.up, this.lnc))
        nbdownlnc[[sp]][paste(tiss,age1,age2,sep="_")]=length(intersect(this.down, this.lnc))
      }
    }
  }

  chicken.tissue.order=setdiff(tissue.order, "Testis")
  chicken.age.order=c("EarlyEmbryo", "LateEmbryo")

  prepare=FALSE
}

###################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure7.pdf", sep=""), width=6.85, height=5)

###################################################################################

m=matrix(rep(NA, 11*9), nrow=11)

for(i in 1:5){
  m[i,]=c(rep(1,3), rep(3, 2), rep(5, 2), rep(7,2))
}

for(i in 6:10){
  m[i,]=c(rep(2,3), rep(4, 2), rep(6, 2), rep(8,2))
}

for(i in 11){
  m[i,]=c(rep(9,9))
}

layout(m)

###################################################################################

## differences between consecutive stages, up-regulated genes

xpos.sp=c(1, 2.5, 4)
names(xpos.sp)=c("Mouse", "Rat", "Chicken")

xadd.tissue=c(-0.45, -0.15, 0.15, 0.45)
names(xadd.tissue)=tissue.order

tinyx=0.08

for(age in c("EarlyEmbryo_LateEmbryo", "LateEmbryo_Newborn", "Newborn_Adult", "Adult_Aged")){
  
  if(age=="EarlyEmbryo_LateEmbryo"){
  
    xlim=c(0.5,5)
  } else{
    xlim=c(0.5,3.75)
  }

  max.up=max(c(unlist(lapply(nbuppc, function(x) max(as.numeric(x[grep(age, names(x), value=T)])))), unlist(lapply(nbuplnc, function(x) max(as.numeric(x[grep(age, names(x), value=T)]))))))
  max.down=max(c(unlist(lapply(nbdownpc, function(x) max(as.numeric(x[grep(age, names(x), value=T)])))), unlist(lapply(nbdownlnc, function(x) max(as.numeric(x[grep(age, names(x), value=T)]))))))

  ylim.up=c(0, max.up*1.1)
  ylim.down=c(-max.down*1.1, 0)

  ## up-regulated genes
  if(age=="EarlyEmbryo_LateEmbryo"){
    par(mar=c(0.25, 3.1, 2.1, 0.25))
  } else{
    par(mar=c(0.25, 1.1, 2.1, 0.25))
  }
  plot(1, type="n", xlab="",  ylab="", axes=F, xlim=xlim, ylim=ylim.up)
  
  for(sp in c("Mouse", "Rat", "Chicken")){
    for(tiss in c("Brain","Kidney", "Liver", "Testis")){
      this.sample=paste(tiss, age, sep="_")
      
      if(this.sample%in%names(nbuppc[[sp]])){
        this.xpos=xpos.sp[sp]+xadd.tissue[tiss]
        this.up.pc=nbuppc[[sp]][[this.sample]]
        this.up.lnc=nbuplnc[[sp]][[this.sample]]
        
        ## points(this.xpos, this.up.pc, pch=22, col=col.tissues[tiss], bg=col.tissues[tiss])
        ## points(this.xpos, this.up.lnc, pch=22, col=col.tissues[tiss], bg="white")


        rect(this.xpos, 0, this.xpos+tinyx, this.up.pc, col=col.tissues[tiss], border=col.tissues[tiss])
        rect(this.xpos+1.5*tinyx, 0, this.xpos+2.5*tinyx, this.up.lnc, col="white", border=col.tissues[tiss])
      }
    }
  }

  axis(side=2, cex.axis=0.75, mgp=c(3, 0.75,0))

  mtext("mouse", side=1, line=-0.25, at=xpos.sp["Mouse"], xpd=NA, cex=0.6)
  mtext("rat", side=1, line=-0.25, at=xpos.sp["Rat"], xpd=NA, cex=0.6)
  
  if(age=="EarlyEmbryo_LateEmbryo"){
    mtext("chicken", side=1, line=-0.25, at=xpos.sp["Chicken"], xpd=NA, cex=0.6)

    mtext("nb. up-regulated genes", side=2, cex=0.6, line=2)
  }

  ## down-regulated genes

  ## up-regulated genes
  if(age=="EarlyEmbryo_LateEmbryo"){
    par(mar=c(2.1, 3.1, 0.25, 0.25))
 } else{
   par(mar=c(2.1,  1.1, 0.25, 0.25))
 }
  
  plot(1, type="n", xlab="",  ylab="", axes=F, xlim=xlim, ylim=ylim.down)

   for(sp in c("Mouse", "Rat", "Chicken")){
    for(tiss in c("Brain","Kidney", "Liver", "Testis")){
      this.sample=paste(tiss, age, sep="_")
      
      if(this.sample%in%names(nbuppc[[sp]])){
        this.xpos=xpos.sp[sp]+xadd.tissue[tiss]
        this.down.pc=nbdownpc[[sp]][[this.sample]]
        this.down.lnc=nbdownlnc[[sp]][[this.sample]]
        
        ## points(this.xpos, -this.down.pc, pch=22, col=col.tissues[tiss], bg=col.tissues[tiss])
        ## points(this.xpos, -this.down.lnc, pch=22, col=col.tissues[tiss], bg="white")

        rect(this.xpos, 0, this.xpos+tinyx,   -this.down.pc,  col=col.tissues[tiss], border=col.tissues[tiss])
        rect(this.xpos+1.5*tinyx, 0, this.xpos+2.5*tinyx, -this.down.lnc,  col="white", border=col.tissues[tiss])
      }
    }
  }

  if(age=="EarlyEmbryo_LateEmbryo"){
    xpos.text=2.5
    ypos.text=ylim.down[1]
    
    text("mid-stage\nembryo", x=xpos.text[1]-0.85, y=ypos.text, cex=0.85, adj=0.5, xpd=NA)
    text("late\nembryo", x=xpos.text[1]+1.35, y=ypos.text, cex=0.85, adj=0.5, xpd=NA)
    arrows(xpos.text[1]+0.2, ypos.text, xpos.text[1]+0.6, ypos.text, length=0.03)

    mtext("nb. down-regulated genes", side=2, cex=0.6, line=2)
  }

  
  if(age=="LateEmbryo_Newborn"){
    xpos.text=1.65
    ypos.text=ylim.down[1]
    
    text("late\nembryo", x=xpos.text[1]-0.6, y=ypos.text, cex=0.85, adj=0.5, xpd=NA)
    text("newborn", x=xpos.text[1]+1.1, y=ypos.text, cex=0.85, adj=0.5, xpd=NA)
    arrows(xpos.text[1]+0.1, ypos.text, xpos.text[1]+0.5, ypos.text, length=0.03)
  }


  
  if(age=="Newborn_Adult"){
    xpos.text=1.65
    ypos.text=ylim.down[1]
    
    text("newborn", x=xpos.text[1]-0.6, y=ypos.text, cex=0.85, adj=0.5, xpd=NA)
    text("young\nadult", x=xpos.text[1]+1.15, y=ypos.text, cex=0.85, adj=0.5, xpd=NA)
    arrows(xpos.text[1]+0.1, ypos.text, xpos.text[1]+0.5, ypos.text, length=0.03)
  }


  
  if(age=="Adult_Aged"){
    xpos.text=1.65
    ypos.text=ylim.down[1]
    
    text("young\nadult", x=xpos.text[1]-0.6, y=ypos.text, cex=0.85, adj=0.5, xpd=NA)
    text("aged\nadult", x=xpos.text[1]+1.1, y=ypos.text, cex=0.85, adj=0.5, xpd=NA)
    arrows(xpos.text[1]+0.1, ypos.text, xpos.text[1]+0.5, ypos.text, length=0.03)
  }
  
  yax=pretty(ylim.down)
  yax=yax[which(yax>=ylim.down[1] & yax<=ylim.down[2])]
  
  axis(side=2, cex.axis=0.75, mgp=c(3, 0.75,0), at=yax, labels=-yax)

}

###################################################################################


## legend

par(mar=c(1.1,1.1,0.1,0.1))
     
plot(1, type="n", xlab="",ylab="", axes=F)

legend("topleft", legend=c("protein-coding", "lncRNA"), pch=22, col="black", pt.bg=c("black", "white"), inset=c(0.055, -0.25), bty="n", cex=0.9, pt.cex=1.4, xpd=NA, horiz=T)

legend("topleft", legend=shortname.tiss, fill=col.tissues, border=col.tissues, bty="n", inset=c(0.05, 0.2), cex=0.9, horiz=T, xpd=NA)

###################################################################################

dev.off()

###################################################################################
