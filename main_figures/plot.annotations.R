######################################################################################
######################################################################################

plot.annot <- function(annot, start, end, plot.axis=FALSE, strands=c("1", "-1"), col.fwd="gray30", col.rev="gray50", biotypes=c("protein_coding"), plot.names=T){

  annot=annot[which(annot$biotype%in%biotypes),]
  annot=annot[order(annot$start),]
  annot=annot[order(annot$geneid),]
  annot=annot[order(annot$strand),]
  
  ## if both strands
  
  ylim=c(0.05,1.2)
  
  if(length(strands)==1 & (strands[1]=="1" | strands[1]=="+")){
    ylim=c(0.6,1.1)
  }
  
  if(length(strands)==1 &  (strands[1]=="-1" | strands[1]=="-")){
    ylim=c(0.2,0.6)
  }

  annot=annot[which(annot$strand%in%strands),]
  

  xlim=c(start, end)
  plot(1,type="n",xlab="",ylab="",axes=F,xlim=xlim,ylim=ylim,xaxs="i")

  genes=unique(annot$geneid)
  ## print(genes)

  height=0.04
  small=0.035

  alltx.fwd=unique(annot$txid[which(annot$strand%in%c("1", "+"))])
  alltx.rev=unique(annot$txid[which(annot$strand%in%c("-1", "-"))])

  propfwd=0
  if(length(alltx.fwd)>0){
    propfwd=length(alltx.fwd)/(length(alltx.fwd)+length(alltx.rev))
    ypos.fwd=seq(from=proprev+0.05, to=1.05, length=length(alltx.fwd))
    names(ypos.fwd)=alltx.fwd
  }
  
  proprev=0

  if(length(alltx.rev)>0){
    proprev=length(alltx.rev)/(length(alltx.fwd)+length(alltx.rev))
    ypos.rev=seq(from=0.05, to=proprev, length=length(alltx.rev))
    names(ypos.rev)=alltx.rev
  }
  
  col=c(col.fwd,col.rev,col.fwd, col.rev)
  names(col)=c("1", "-1","+", "-")
 
  if(plot.axis){
    xaxis=pretty(c(start,end))
    xlabels=paste(round(xaxis/1000, digits=0),"kb")
    axis(side=1, at=xaxis, labels=xlabels, cex.axis=0.85)
  }
  arrow.width=diff(xlim)/50
  
  for(g in genes){
    print(g)

    tx=unique(annot$txid[which(annot$geneid==g)])

    for(t in tx){
      
      this.coords=annot[which(annot$geneid==g & annot$txid==t),]
      this.coords=this.coords[order(this.coords$start),]
      
      nbex=dim(this.coords)[1]
      
      strand=as.character(this.coords$strand[1])
      name=this.coords$genename[1]
      
      if(strand=="1" | strand=="+"){
        this.ypos=ypos.fwd[t]
      }

      if(strand=="-1" | strand=="-"){
        this.ypos=ypos.rev[t]
      }
      
      this.col=col[strand]
      
      minx=min(this.coords$start)
      maxx=max(this.coords$end)
      
      for(i in 1:nbex){
        thisstart=this.coords[i,"start"]
        thisend=this.coords[i,"end"]
        
        rect(thisstart,this.ypos,thisend,this.ypos+height,col=this.col,border=this.col)
        thisy=this.ypos+height
        ymid=thisy+height/2
        
        if(i>1){
          prevstart=this.coords[i-1,"start"]
          prevend=this.coords[i-1,"end"]
          
          midpoint=mean(c(prevend,thisstart))
          segments(prevend,thisy,midpoint,ymid,col=this.col)
          segments(midpoint,ymid,thisstart,thisy,col=this.col)
        }
      }
      
      if(strand=="1" | strand=="+"){
        tinyx=(end-start)/1000

        if(plot.names){
          ## text(x=minx-tinyx,y=this.ypos+4*small,labels=name,cex=0.9, adj=0, font=3, xpd=NA) ## names of the genes
          text(x=minx-5*tinyx,y=this.ypos+height/2,labels=name,cex=0.9, adj=1, font=3, xpd=NA) ## names of the genes
        }
        
        segments(minx,this.ypos+small, minx, this.ypos+height+small*1.3, col=this.col)
        arrows(minx, this.ypos+height+small*1.3, minx+arrow.width, this.ypos+height+small*1.3, length=0.05, col=this.col)
      }
      
      if(strand=="-1" | strand=="-"){
        tinyx=(end-start)/1000
        
        if(plot.names){
       ##  text(x=maxx-tinyx,y=this.ypos+4*small,labels=name,cex=0.9, adj=1, font=3, xpd=NA) ## names of the genes
          text(x=maxx+5*tinyx,y=this.ypos+height/2,labels=name,cex=0.9, adj=0, font=3, xpd=NA) ## names of the genes
        }
        
        segments(maxx,this.ypos+small, maxx, this.ypos+height+small*1.3, col=this.col)
        arrows(maxx, this.ypos+height+small*1.3, maxx-arrow.width, this.ypos+height+small*1.3, length=0.05, col=this.col)
      }
    }
  }
}

######################################################################################
######################################################################################

plot.annot.flat <- function(annot, start, end, plot.axis=FALSE, strands=c("1", "-1"), col.fwd="gray30", col.rev="gray50", biotypes=c("protein_coding"), ylim=NA, ypos=NA){

  annot=annot[which(annot$biotype%in%biotypes),]
  annot=annot[order(annot$start),]
  annot=annot[order(annot$geneid),]
  annot=annot[order(annot$strand),]
  
  annot=annot[which(annot$strand%in%strands),]
  
  ## if both strands

  if(any(is.na(ylim))){
    ylim=c(0.05,1.2)
    
    if(length(strands)==1 & (strands[1]=="1" | strands[1]=="+")){
      ylim=c(0.6,1.1)
    }
    
    if(length(strands)==1 &  (strands[1]=="-1" | strands[1]=="-")){
      ylim=c(0.2,0.6)
    }
  }
  

  xlim=c(start, end)
  plot(1,type="n",xlab="",ylab="",axes=F,xlim=xlim,ylim=ylim,xaxs="i")

  genes=unique(annot$geneid)
  ## print(genes)

  height=0.04
  small=0.05

  allgene.fwd=unique(annot$geneid[which(annot$strand%in%c("1", "+"))])
  allgene.rev=unique(annot$geneid[which(annot$strand%in%c("-1", "-"))])

  propfwd=length(allgene.fwd)/(length(allgene.fwd)+length(allgene.rev))
  proprev=length(allgene.rev)/(length(allgene.fwd)+length(allgene.rev))

  if(any(is.na(ypos))){
    ypos.fwd=seq(from=proprev+0.05, to=1.05, length=length(allgene.fwd))
    ypos.rev=seq(from=0.05, to=proprev, length=length(allgene.rev))
  } else{
    ypos.fwd=rep(ypos, length(allgene.fwd))
    ypos.rev=rep(ypos, length(allgene.rev))
  }
  
  names(ypos.fwd)=allgene.fwd
  names(ypos.rev)=allgene.rev

  col=c(col.fwd,col.rev,col.fwd, col.rev)
  names(col)=c("1", "-1","+", "-")
 
  if(plot.axis){
    xaxis=pretty(c(start,end))
    xlabels=paste(round(xaxis/1000, digits=0),"kb")
    axis(side=1, at=xaxis, labels=xlabels, cex.axis=0.85)
  }
  arrow.width=diff(xlim)/50
  
  for(g in genes){
    print(g)
    
    this.coords=annot[which(annot$geneid==g),]
    this.coords=this.coords[order(this.coords$start),]
    
    nbex=dim(this.coords)[1]

    strand=as.character(this.coords$strand[1])
    
    if(strand=="1" | strand=="+"){
      this.ypos=ypos.fwd[g]
    }
    
    if(strand=="-1" | strand=="-"){
      this.ypos=ypos.rev[g]
    }
    
    name=this.coords$genename[1]
    
    this.col=col[strand]
    
    minx=min(this.coords$start)
    maxx=max(this.coords$end)

    if(strand=="1" | strand=="+"){
      tinyx=(end-start)/1000
      #text(x=minx-tinyx,y=this.ypos+5*small,labels=name,cex=0.9, adj=0, font=3, xpd=NA) ## names of the genes
      
      text(x=minx-5*tinyx,y=this.ypos+height/2,labels=name,cex=0.9, adj=1, font=3, xpd=NA) ## names of the genes
      
      segments(minx,this.ypos+small, minx, this.ypos+height+small*1.3, col=this.col)
      arrows(minx, this.ypos+height+small*1.3, minx+arrow.width, this.ypos+height+small*1.3, length=0.05, col=this.col, xpd=NA)
    }
    
    if(strand=="-1" | strand=="-"){
      tinyx=(end-start)/1000
     ##  text(x=maxx-tinyx,y=this.ypos+5*small,labels=name,cex=0.9, adj=1, font=3, xpd=NA) ## names of the genes
      
       text(x=maxx+5*tinyx,y=this.ypos+height/2,labels=name,cex=0.9, adj=0, font=3, xpd=NA) ## names of the genes
      
      segments(maxx,this.ypos+small, maxx, this.ypos+height+small*1.3, col=this.col)
      arrows(maxx, this.ypos+height+small*1.3, maxx-arrow.width, this.ypos+height+small*1.3, length=0.05, col=this.col, xpd=NA)
    }

    
    for(i in 1:nbex){
      thisstart=this.coords[i,"start"]
      thisend=this.coords[i,"end"]

      rect(thisstart,this.ypos,thisend,this.ypos+height,col=this.col,border=this.col)
      thisy=this.ypos+height
      ymid=thisy+height/2
      
      if(i>1){
        prevstart=this.coords[i-1,"start"]
        prevend=this.coords[i-1,"end"]
        
        segments(prevend,this.ypos+height/2,thisstart,this.ypos+height/2,col=this.col)
      }
    }
  }
}

######################################################################################
######################################################################################
