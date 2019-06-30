##############################################################################

objects=ls()

if(!("pathFigures"%in%objects)){
  source("parameters.R")
  load=TRUE
  prepare=TRUE
}

##############################################################################

if(load==TRUE){

  load("RData/data.expression.ortho.RData")

  avgexp.mr=avgexp.mr[,-c(1,2)]

  samples.mr=colnames(avgexp.mr)
  species.mr=unlist(lapply(colnames(avgexp.mr), function(x) unlist(strsplit(x,split="_"))[1]))
  tissue.mr=unlist(lapply(colnames(avgexp.mr), function(x) unlist(strsplit(x,split="_"))[2]))
  age.mr=unlist(lapply(colnames(avgexp.mr), function(x) unlist(strsplit(x,split="_"))[3]))
  tissage.mr=paste(tissue.mr, age.mr, sep="_")

  ## organ/stage markers

  markers=read.table(paste(pathTables, "SupplementaryTable4.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

  go.enrichment=list()

  for(tiss in tissue.order){
    for(age in age4.order){
      this.path=paste(pathDatasets, "SupplementaryDataset3/GOEnrichment_BiologicalProcess_OrganStageMarkers_MouseRat_",tiss,"_",age,".txt", sep="")

      if(file.exists(this.path)){
        this.go=read.table(this.path, h=T, sep="\t", quote="\"",stringsAsFactors=F)

        go.enrichment[[paste(tiss, age, sep="_")]]=this.go
      }
    }
  }
  
  load=FALSE
}

##############################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

indexfigures=c(2,3,4,5)
names(indexfigures)=tissue.order

heights=c(8.25,8.25,8.25,6.25)
names(heights)=tissue.order

nbplots=c(4,4,4,3)
names(nbplots)=tissue.order

labels=c("A", "B", "C", "D", "E","F","G","H")

nbmaxgo=7

maxchar=50

###################################################################################

for(tiss in tissue.order){

  ###################################################################################
  
  pdf(file=paste(pathFigures, "SupplementaryFigure",indexfigures[tiss],".pdf", sep=""), height=heights[tiss], width=6.85)

  ###################################################################################
  
  m=matrix(rep(NA, 2*nbplots[tiss]), nrow=2)

  m[1,]=seq(from=1, to=(2*nbplots[tiss]-1), by=2)
  m[2,]=seq(from=2, to=2*nbplots[tiss], by=2)

  m=t(m)

  layout(m)

##############################################################################
  
  xpos=1:20
  names(xpos)=kronecker(tissue.order, age.order, paste, sep="_")
  
  xsp=c(-0.2, 0.2)
  names(xsp)=c("Mouse", "Rat")

  density=c(200,20)
  names(density)=names(xsp)
  
  xlim=c(0.5,20.5)
  
  width=0.1

  indexplot=0
  
  for(age in age4.order){
    nbgenes=length(which(markers$MaxSample.Mouse==paste(tiss, age, sep="_")))
    
    if(nbgenes>0){
      indexplot=indexplot+1

      label=labels[indexplot]
      
      gene=paste(markers[which(markers$MaxSample.Mouse==paste(tiss, age, sep="_"))[1],c("ID.Mouse", "ID.Rat")], collapse="_")
      name=markers[which(markers$MaxSample.Mouse==paste(tiss, age, sep="_"))[1],"GeneName"]
      this.exp=as.numeric(avgexp.mr[gene,])
      names(this.exp)=samples.mr
      
      ylim=c(0, max(this.exp))
      
      par(mar=c(3.25,3.1,1.5,0.5))
            
      plot(1, type="n", xlim=xlim, ylim=ylim, axes=F, xlab="", ylab="")
      rect(xpos[tissage.mr]+xsp[species.mr]-width,0,  xpos[tissage.mr]+xsp[species.mr]+width, this.exp, col=col.tissues[tissue.mr], border=col.tissues[tissue.mr], density=density[species.mr])
      
      axis(side=1, at=xpos, labels=rep("", 20),mgp=c(3,0.5,0))
      axis(side=2, mgp=c(3,0.5,0))

      mtext("expression level (TPM)", side=2, cex=0.65, line=2)
      mtext(rep(1:5, 4), at=xpos, cex=0.55, side=1, line=0.5)
      
      mtext(name, side=3, line=0.1, cex=0.65, font=3)
      mtext(label, side=3, line=0.45, cex=0.85, font=2, at=-2.85)
    
      mtext("developmental stage", side=1, cex=0.65, line=1.5)
      
      if(label=="A"){
        if(tiss=="Testis"){
          legend("topleft", legend=shortname.tiss, fill=col.tissues, border=col.tissues, bty="n", inset=c(0.05,-0.01), cex=0.95, horiz=F, xpd=NA)
        } else{
          legend("topright", legend=shortname.tiss, fill=col.tissues, border=col.tissues, bty="n", inset=c(0.05,-0.01), cex=0.95, horiz=F, xpd=NA)
        }
      }

      if(label=="C"){
        if(tiss=="Testis"){
          legend("topleft", legend=c("mouse", "rat"), fill="black", density=density, bty="n", inset=c(0.05,-0.01), cex=0.95, horiz=F, xpd=NA)
        } else{
          legend("topright", legend=c("mouse", "rat"), fill="black", density=density, bty="n", inset=c(0.05,-0.01), cex=0.95, horiz=F, xpd=NA)
        }
      }
      
            
      ## gene ontology

      if(paste(tiss, age, sep="_")%in%names(go.enrichment)){
        indexplot=indexplot+1

        label=labels[indexplot]
        
        this.go=go.enrichment[[paste(tiss, age, sep="_")]]
       
        if(dim(this.go)[1]>nbmaxgo){
          this.go=this.go[1:nbmaxgo,]
        }

        nbchar=nchar(this.go$GOName)
        this.go$GOName[which(nbchar>maxchar)]=unlist(lapply(this.go$GOName[which(nbchar>maxchar)], function(x) {y=unlist(strsplit(x, split=" ")); z=paste(y[1:floor(length(y)/2)], collapse=" "); t=paste(y[(floor(length(y)/2)+1):length(y)], collapse=" "); u=paste(z, t, sep="\n"); return(u); }))
        
        ylim.go=c(0.5,nbmaxgo+0.5)
        this.go$FDR[which(this.go$FDR<1e-10)]=1e-10
        xlim.go=c(0, 10)
        smallx=0.5
        
        par(mar=c(3.25,1.1,1.5,20))
        
        plot(1, type="n",  axes=F, xlab="", ylab="", xlim=xlim.go, ylim=ylim.go, xaxs="i", yaxs="i")

        axis(side=1, mgp=c(3,0.5,0), cex.axis=0.85, at=c(0,5,10))
      
        mtext("-log10(FDR)", cex=0.65, side=1, line=1.5)
        
        nbgo=dim(this.go)[1]
        ypos=(nbmaxgo:1)[1:nbgo]

        axis(side=2, at=nbmaxgo:1, labels=rep("", nbmaxgo))
        
        ywidth=0.15
        
        rect(0, ypos-ywidth, -log10(this.go$FDR),ypos+ywidth, col="gray60")
        for(i in 1:dim(this.go)[1]){
          text(x=-log10(this.go$FDR[i])+smallx, y=ypos[i], labels=this.go$GOName[i],adj=c(0,0.5), xpd=NA)
        }
        
        mtext(label, side=3, line=0.45, cex=0.85, font=2, at=-0.5)
        
      } else{
        plot(1, type="n",  axes=F, xlab="", ylab="")
      }
    }
  }
  
  ##############################################################################
  
  dev.off()
}

##############################################################################

