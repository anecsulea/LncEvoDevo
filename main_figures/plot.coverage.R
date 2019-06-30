######################################################################################
######################################################################################

plot.coverage <- function(cov,xstart,xend, plot.axis=TRUE, recompute.ylim=TRUE, ylim=NA, col="black", mtext.xaxis=""){
  cov=cov[which((cov$V2>=xstart & cov$V2<=xend) | (cov$V3>=xstart & cov$V3<=xend)),]
  nbcov=dim(cov)[1]

 
  if(nbcov>0){
    if(recompute.ylim & length(cov)>0){
      ylim=c(0,max(cov$V4))
    }
    
    if(recompute.ylim & length(cov)==0){
      ylim=c(0,1)
    }
    
    plot(1,type="n",xlab="",ylab="",axes=F,xlim=c(xstart,xend),ylim=ylim,xaxs="i")
    
    yaxis=pretty(ylim,n=2)
    yaxis=yaxis[which(yaxis<=ylim[2])]
    
    axis(side=2, cex.axis=0.5, mgp=c(3,0.5,0), at=yaxis, labels=rep("", length(yaxis)))
    mtext(yaxis, at=yaxis, side=2, cex=0.6, line=0.5)
    
    ## for the x axis 
    if(plot.axis){
     xaxis=pretty(c(xstart,xend))
     xlabels=paste(round(xaxis/1000, digits=0),"kb")
     axis(side=1, at=xaxis, labels=xlabels,cex.axis=0.85, mgp=c(3,0.35,0))
     ## add mtext if needed
     
     mtext(mtext.xaxis, side=1, at=xstart, line=0.37, cex=0.6)
    } else{
      xaxis=pretty(c(xstart,xend))
      xlabels=rep("", length(xaxis))
      axis(side=1, at=xaxis, labels=xlabels,cex.axis=0.85, mgp=c(3,0.35,0))
    }
      
    if(length(cov)>0){
      rect(cov$V2+1,0,cov$V3,cov$V4,col=col,border=col)
    }
    
  }
  else{
     plot(1,type="n",xlab="",ylab="",axes=F,xlim=c(xstart,xend),ylim=ylim,xaxs="i")
     axis(side=2)
     xaxis=pretty(c(start,end))
     xlabels=paste(round(xaxis/1000, digits=0),"kb")
     axis(side=1, at=xaxis, labels=xlabels)
  }
}

########################################################################
########################################################################
