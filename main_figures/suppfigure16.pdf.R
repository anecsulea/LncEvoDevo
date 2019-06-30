###################################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE
}

###################################################################################

if(load==TRUE){
   for(sp in c("Mouse", "Rat")){
     load(paste("RData/data.annotations.", sp, ".RData", sep=""))
     assign(paste("pc", tolower(sp), sep="."), pc)
     assign(paste("lnc", tolower(sp), sep="."), lnc)
     assign(paste("allinfo", tolower(sp), sep="."), allinfo)
     
     load(paste("RData/data.expression.",sp,".RData", sep=""))
     assign(paste("normtpm", tolower(sp), sep="."), normtpm)

     samples=colnames(normtpm)
     tissue=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))
     age=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[2]))
     age=unlist(lapply(age, function(x) substr(x, 1, nchar(x)-1)))

     assign(paste("samples", tolower(sp), sep="."), samples)
     assign(paste("tissue", tolower(sp), sep="."), tissue)
     assign(paste("age", tolower(sp), sep="."), age)
     
     rm(list=c("pc","lnc", "allinfo", "rawtpm", "normtpm", "kcounts", "readcounts", "downsampled", "samples", "tissue", "age"))
   }

   
   load("RData/data.comparison.specific.ortho.RData")
   
   load=FALSE
}

###################################################################################

if(prepare==TRUE){
  
  for(sp in c("Mouse", "Rat")){
    this.samples=get(paste("samples", tolower(sp), sep="."))
    this.tissue=get(paste("tissue", tolower(sp), sep="."))
    this.age=get(paste("age", tolower(sp), sep="."))
    this.normtpm=get(paste("normtpm", tolower(sp), sep="."))
    this.pc=get(paste("pc", tolower(sp), sep="."))
    this.lnc=get(paste("lnc", tolower(sp), sep="."))

    maxexp=log2(apply(this.normtpm, 1, function(x) max(x, na.rm=T))+1)
    maxtissue=apply(this.normtpm, 1, function(x) this.tissue[which.max(x)])
    maxage=apply(this.normtpm, 1, function(x) this.age[which.max(x)])

    maxtissue[which(maxexp==0)]=NA
    maxage[which(maxexp==0)]=NA

    names(maxtissue)=rownames(this.normtpm)
    names(maxage)=rownames(this.normtpm)
    names(maxexp)=rownames(this.normtpm)
        
    assign(paste("maxage", tolower(sp), sep="."), maxage)
    assign(paste("maxtissue", tolower(sp), sep="."), maxtissue)
    assign(paste("maxexp", tolower(sp), sep="."), maxexp)

  }

  ## tissue/stage specificity

  for(sp in c("Mouse", "Rat", "Chicken")){
    load(paste("RData/data.tpm.stats.",sp,".RData", sep=""))

    assign(paste("stats", tolower(sp), sep="."), stats) 
  }
  
  prepare=FALSE
}

###################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################


pdf(file=paste(pathFigures, "SupplementaryFigure16.pdf", sep=""), width=6.85, height=2.5)

m=matrix(rep(NA, 10*10), nrow=10)

for(i in 1:9){
  m[i,]=c(rep(1:5, each=2))
}


for(i in c(10)){
  m[i,]=c(rep(6,10))
}

layout(m)

###################################################################################

## proportion of overlap with Encode enhancers

prop.ov.specific.new=100*length(which(allinfo.mouse[intersect(mouse.specific, mouse.new), "OverlapEncodeEnhancer"]=="Yes"))/length(which(allinfo.mouse[intersect(mouse.specific, mouse.new), "OverlapEncodeEnhancer"]%in%c("Yes", "No")))

prop.ov.ortho.new=100*length(which(allinfo.mouse[intersect(mouse.ortho, mouse.new), "OverlapEncodeEnhancer"]=="Yes"))/length(which(allinfo.mouse[intersect(mouse.ortho, mouse.new), "OverlapEncodeEnhancer"]%in%c("Yes", "No")))

prop.ov.specific.known=100*length(which(allinfo.mouse[intersect(mouse.specific, mouse.known), "OverlapEncodeEnhancer"]=="Yes"))/length(which(allinfo.mouse[intersect(mouse.specific, mouse.known), "OverlapEncodeEnhancer"]%in%c("Yes", "No")))

prop.ov.ortho.known=100*length(which(allinfo.mouse[intersect(mouse.ortho, mouse.known), "OverlapEncodeEnhancer"]=="Yes"))/length(which(allinfo.mouse[intersect(mouse.ortho, mouse.known), "OverlapEncodeEnhancer"]%in%c("Yes", "No")))

allprop=c(prop.ov.ortho.new, prop.ov.specific.new, prop.ov.ortho.known, prop.ov.specific.known)

par(mar=c(3.5,3.1,2.1,0.5))
xpos=c(1:2, 3.5,4.5)
plot(1, type="n", xlim=c(0.5,max(xpos)+0.5), ylim=c(0,55), axes=F, xlab="", ylab="")

width=0.35
rect(xpos-width, 0,xpos+width, allprop, col=c("steelblue4", "aliceblue"))

axis(side=2, cex.axis=0.85, mgp=c(3,0.75,0))
axis(side=1, cex.axis=0.75, mgp=c(3,0.5,0), labels=rep("",4), at=xpos)

mtext("% overlap with enhancers", side=2, line=1.95, cex=0.65)

mtext(rep(c("1-1", "sp"),2), side=1, at=xpos, line=0.5, cex=0.6)
mtext(c("new", "Ensembl"), at=c(1.5,4),side=1, line=1.5, cex=0.6) 

mtext("A", side=3, line=1, at=-1.3, font=2, cex=0.9)


#####################################################################################

## proportion spliced genes

prop.spliced.mouse.specific=100*length(which(allinfo.mouse[mouse.specific, "NbExons"]>=2))/length(mouse.specific)
prop.spliced.mouse.ortho=100*length(which(allinfo.mouse[mouse.ortho, "NbExons"]>=2))/length(mouse.ortho)

prop.spliced.rat.specific=100*length(which(allinfo.rat[rat.specific, "NbExons"]>=2))/length(rat.specific)
prop.spliced.rat.ortho=100*length(which(allinfo.rat[rat.ortho, "NbExons"]>=2))/length(rat.ortho)
  
allprop=c(prop.spliced.mouse.ortho, prop.spliced.mouse.specific, prop.spliced.rat.ortho, prop.spliced.rat.specific)

par(mar=c(3.5,3.1,2.1,1.1))
xpos=c(1:2, 3.5,4.5)
plot(1, type="n", xlim=c(0.5,max(xpos)+0.5), ylim=c(0,65), axes=F, xlab="", ylab="")

width=0.35
rect(xpos-width, 0,xpos+width, allprop, col=c("steelblue4", "aliceblue"))

axis(side=2, cex.axis=0.85, mgp=c(3,0.75,0))
axis(side=1, cex.axis=0.75, mgp=c(3,0.5,0), labels=rep("",4), at=xpos)

mtext("% multi-exonic loci", side=2, line=1.95, cex=0.65)
mtext(rep(c("1-1", "sp"),2), side=1, at=xpos, line=0.5, cex=0.6)
mtext(c("mouse", "rat"), at=c(1.5,4),side=1, line=1.5, cex=0.6) 

mtext("B", side=3, line=1, at=-1.3, font=2, cex=0.9)

## #####################################################################################


prop.bidir.mouse.specific=100*length(which(!is.na(allinfo.mouse[mouse.specific, "BidirectionalPromoter"])))/length(mouse.specific)
prop.bidir.mouse.ortho=100*length(which(!is.na(allinfo.mouse[mouse.ortho, "BidirectionalPromoter"])))/length(mouse.ortho)

prop.bidir.rat.specific=100*length(which(!is.na(allinfo.rat[rat.specific, "BidirectionalPromoter"])))/length(rat.specific)
prop.bidir.rat.ortho=100*length(which(!is.na(allinfo.rat[rat.ortho, "BidirectionalPromoter"])))/length(rat.ortho)
  
allprop=c(prop.bidir.mouse.ortho, prop.bidir.mouse.specific, prop.bidir.rat.ortho, prop.bidir.rat.specific)

par(mar=c(3.5,3.1,2.1,1.1))
xpos=c(1:2, 3.5,4.5)
plot(1, type="n", xlim=c(0.5,max(xpos)+0.5), ylim=c(0,65), axes=F, xlab="", ylab="")

width=0.35
rect(xpos-width, 0,xpos+width, allprop, col=c("steelblue4", "aliceblue"))

axis(side=2, cex.axis=0.85, mgp=c(3,0.75,0))
axis(side=1, cex.axis=0.75, mgp=c(3,0.5,0), labels=rep("",4), at=xpos)

mtext("% bidirectional promoters", side=2, line=1.95, cex=0.65)
mtext(rep(c("1-1", "sp"),2), side=1, at=xpos, line=0.5, cex=0.6)
mtext(c("mouse", "rat"), at=c(1.5,4),side=1, line=1.5, cex=0.6) 

mtext("C", side=3, line=1, at=-1.3, font=2, cex=0.9)

#######################################################################################


## table for tissue 
tab.mouse.specific=as.numeric(table(factor(maxtissue.mouse[mouse.specific], levels=tissue.order)))
tab.mouse.ortho=as.numeric(table(factor(maxtissue.mouse[mouse.ortho], levels=tissue.order)))

tab.rat.specific=as.numeric(table(factor(maxtissue.rat[rat.specific], levels=tissue.order)))
tab.rat.ortho=as.numeric(table(factor(maxtissue.rat[rat.ortho], levels=tissue.order)))


m=matrix(c(tab.mouse.ortho,tab.mouse.specific, tab.rat.ortho,tab.rat.specific), nrow=4, byrow=T)
m=m/apply(m,1,sum)

xpos=c(1, 1.35, 2, 2.35)

width=0.15

par(mar=c(3.5,3,2.1,0.5))

plot(1, xlim=c(0.8,2.55), ylim=c(0,100), axes=F, xlab="", ylab="", type="n")

for(i in 1:4){
  ypos=100*cumsum(m[i,])
  rect(xpos[i]-width, c(0, ypos[-length(ypos)]), xpos[i]+width, ypos, col=col.tissues, border=NA)
}

axis(side=2, cex.axis=0.8, mgp=c(3,0.75,0))
axis(side=1, at=xpos, labels=rep("",4))

mtext(rep(c("1-1", "sp"),2), at=xpos, line=0.5, side=1, cex=0.6)
mtext(c("mouse", "rat"), at=c(mean(xpos[1:2]), mean(xpos[3:4])),side=1, cex=0.6, line=1.5)

mtext("% genes w. max. exp. in organ", side=2, cex=0.65, line=1.9)

mtext("D", side=3, at=0.05, line=1, font=2, cex=0.9)


##########################################################################

## stage with maximum expression


## table for age 
tab.mouse.specific=as.numeric(table(factor(maxage.mouse[mouse.specific], levels=age.order)))
tab.mouse.ortho=as.numeric(table(factor(maxage.mouse[mouse.ortho], levels=age.order)))

tab.rat.specific=as.numeric(table(factor(maxage.rat[rat.specific], levels=age.order)))
tab.rat.ortho=as.numeric(table(factor(maxage.rat[rat.ortho], levels=age.order)))


m=matrix(c(tab.mouse.ortho,tab.mouse.specific, tab.rat.ortho,tab.rat.specific), nrow=4, byrow=T)
m=m/apply(m,1,sum)

xpos=c(1, 1.35, 2, 2.35)

width=0.15

par(mar=c(3.5,3,2.1,0.5))

plot(1, xlim=c(0.8,2.55), ylim=c(0,100), axes=F, xlab="", ylab="", type="n")

for(i in 1:4){
  ypos=100*cumsum(m[i,])
  rect(xpos[i]-width, c(0, ypos[-length(ypos)]), xpos[i]+width, ypos, col=col.ages, border=NA)
}

axis(side=2, cex.axis=0.8, mgp=c(3,0.75,0))
axis(side=1, at=xpos, labels=rep("",4))

mtext(rep(c("1-1", "sp"),2), at=xpos, line=0.5, side=1, cex=0.6)
mtext(c("mouse", "rat"), at=c(mean(xpos[1:2]), mean(xpos[3:4])),side=1, cex=0.6, line=1.5)

mtext("% genes w. max. exp. in stage", side=2, cex=0.65, line=1.9)

mtext("E", side=3, at=0.05, line=1, font=2, cex=0.9)

##########################################################################
## legend

par(mar=c(0,0,0,0))
plot(1, type="n", xlab="", ylab="", axes=F)
legend("topleft", legend=shortname.tiss, border=col.tissues, fill=col.tissues, horiz=T, bty="n", xpd=NA, inset=c(0.7, -0.1), cex=0.85)
legend("topleft", legend=shortname.age, border=col.ages, fill=col.ages, horiz=T, bty="n", xpd=NA, cex=0.85, inset=c(0.25,0.3))
legend("topleft", legend=c("ortho 1-to-1", "species-specific"), fill=c("steelblue4", "aliceblue"), border=c("black"),horiz=T, bty="n", xpd=NA, cex=0.85, inset=c(0.05,-0.1))

##########################################################################

dev.off()

##########################################################################
