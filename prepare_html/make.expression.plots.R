#################################################################

path="/sps/biometr/necsulea/LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathHTML=paste(path, "html/", sep="")

#################################################################

tissue.order=c("Brain", "Kidney", "Liver", "Testis")
age.order=c("EarlyEmbryo", "LateEmbryo", "Newborn", "Adult", "Aged")

col.tissues=c("dodgerblue4", "seagreen", "indianred4","goldenrod")
names(col.tissues)=c("Brain", "Kidney", "Liver","Testis")

xpos=1:20
names(xpos)=kronecker(tissue.order, age.order, paste, sep="_")

#################################################################

args = commandArgs(trailingOnly=TRUE)

sp=args[1]

#################################################################
  
tpm=read.table(paste(pathDatasets, "SupplementaryDataset2/KallistoNormalizedTPM_",sp,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(tpm)=tpm$GeneID
tpm=tpm[,-1]

samples=colnames(tpm)
tissues=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))
ages=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[2]))
ages=unlist(lapply(ages, function(x) substr(x,1, nchar(x)-1)))

geneinfo=read.table(paste(pathDatasets, "SupplementaryDataset1/GeneInfo_",sp,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(geneinfo)=geneinfo$GeneID

pc=geneinfo$GeneID[which(geneinfo$SelectedProteinCoding=="Yes")]
lnc=geneinfo$GeneID[which(geneinfo$CandidateLncRNA=="Yes")]

for(gene in c(pc, lnc)){
  name=geneinfo[gene, "GeneName"]
  
  if(is.na(name)){
    name=""
  }

  if(gene%in%rownames(tpm)){
    this.tpm=log2(as.numeric(tpm[gene,])+1)
    
    ylim=range(this.tpm)
    xlim=c(1,20)
    
    png(paste(pathHTML, "img/", sp, "_", gene,".png",sep=""), width=800, height=400)
    
    par(mar=c(5.1, 5.1, 2.1, 1.1))
    
    plot(1, xlim=xlim, ylim=ylim, type="n", xlab="", ylab="", axes=F)
    box()
    points(xpos[paste(tissues, ages, sep="_")], this.tpm, col=col.tissues[tissues], pch=20, cex=1.5)
    axis(side=1, at=1:20, labels=rep("",20))
    axis(side=2, cex.axis=1.1)
    mtext(rep(1:5, 4), side=1, at=xpos, line=0.75, cex=1.1)
    mtext("samples", side=1, line=3.5, cex=1.25)
    mtext("expression value (log2-transformed TPM)", side=2, line=3, cex=1.25)
    mtext(paste(sp, gene, name, sep=" "), side=3, cex=1.1)
    
    abline(v=c(5.5,10.5,15.5), lty=2, col="gray40")
    
    mtext(c("brain", "kidney", "liver", "testes"), side=1, at=c(3,8,13,18), line=1.85, cex=1.1)
    
    dev.off()
  }
}

#################################################################


