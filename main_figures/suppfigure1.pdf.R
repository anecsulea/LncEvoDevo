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

    load(paste("RData/data.expression.",sp,".RData", sep=""))
    assign(paste("downsampled", tolower(sp), sep="."), downsampled)
    
    rm(list=c("pc","lnc", "allinfo", "rawtpm", "normtpm", "kcounts", "readcounts", "downsampled"))
  }

  load("RData/data.expression.ortho.RData")
  avgexp.mr=avgexp.mr[,setdiff(colnames(avgexp.mr), c("GeneID", "GeneType"))]

  
  load("RData/data.celltype.markers.RData")
  
  load=F
}

###################################################################################

if(prepare==TRUE){
  
  ## select markers

  markers$CellType[which(markers$CellType=="loop_of_Henle_ascending_limb_epithelial_cell")]="loop_of_Henle_epithelial"
  markers$CellType[which(markers$CellType=="oligodendrocyte_precursor_cell")]="oligodendrocyte_precursor"

  markers$CellType[which(markers$CellType=="proximal_straight_tubule_epithelial_cell")]="proximal_straight_tubule\nepithelial cell"
  markers$CellType[which(markers$CellType=="endothelial_cell_of_hepatic_sinusoid")]="endothelial_cell,_hepatic_sinusoid"
  markers$CellType[which(markers$CellType=="collecting_duct_epithelial_cell")]="collecting_duct\nepithelial_cell"

  
  markers=markers[which(!is.na(markers$RatGeneID)),]
  markers$Organ[which(markers$ImmuneCells=="yes")]="ImmuneCells"
  markers=markers[order(paste(markers$Organ, markers$CellType)),]
  markers=markers[c(which(markers$ImmuneCells=="no"), which(markers$ImmuneCells=="yes")),]

  firstgenes=which(markers$CellType%in%c("oligodendrocyte", "oligodendrocyte_precursor"))
  markers=markers[c(firstgenes, setdiff(1:dim(markers)[1], firstgenes)),]
 
  markers$OrthoID=paste(markers$MouseGeneID, markers$RatGeneID,sep="_")
  markers=markers[which(markers$OrthoID%in%rownames(avgexp.mr)),]
  
  sample.order=kronecker(tissue.order, age.order, paste, sep="_")
  combined.sample.order=kronecker(sample.order, c("Mouse", "Rat"), function(x,y) paste(y,x, sep="_"))
    
  all.tpm=log2(avgexp.mr[markers$OrthoID,intersect(combined.sample.order, colnames(avgexp.mr))]+1)
  
  exp.zscore=(all.tpm-apply(all.tpm, 1, mean, na.rm=T))/apply(all.tpm, 1, sd, na.rm=T)

  exp.zscore$Mouse_Testis_EarlyEmbryo=rep(NA, dim(exp.zscore)[1]) ## we add one NA line

 
  exp.zscore=exp.zscore[,intersect(combined.sample.order, colnames(exp.zscore))]
  
  rownames(exp.zscore)=markers$GeneName

  exp.zscore=exp.zscore[dim(exp.zscore)[1]:1,] ## brain at the top of the list

  exp.zscore=t(exp.zscore)
  
  prepare=F
}

###################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

###################################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure1.pdf", sep=""), width=6.85, height=7.5)

###################################################################################

## layout

m=matrix(rep(NA,22*12), nrow=22)

for(i in c(1:20)){
  m[i,]=c(rep(1,12))
}


for(i in c(21:22)){
  m[i,]=c(rep(2,12))
}

layout(m)

###################################################################################

## heatmaps for markers

par(mar=c(0, 1.5, 2.1, 15.2))
image(exp.zscore, axes=F, col=topo.colors(50))

xpos=seq(from=0, to=1, length=dim(exp.zscore)[1])
ypos=seq(from=0, to=1, length=dim(exp.zscore)[2])
ypos=ypos[length(ypos):1]

this.samples=rownames(exp.zscore)
this.species=unlist(lapply(this.samples, function(x) unlist(strsplit(x, split="_"))[1]))
this.tissues=unlist(lapply(this.samples, function(x) unlist(strsplit(x, split="_"))[2]))
this.ages=unlist(lapply(this.samples, function(x) unlist(strsplit(x, split="_"))[3]))

label.ages=1:5
names(label.ages)=age.order

## labels for tissues

for(tiss in tissue.order){
  mtext(tolower(tiss), side=3, at=mean(xpos[which(this.tissues==tiss)]), line=0.8, cex=0.65, col=col.tissues[tiss], font=2)

  pos.ages=tapply(xpos[which(this.tissues==tiss)], as.factor(this.ages[which(this.tissues==tiss)]), mean)
  names(pos.ages)=levels(as.factor(this.ages[which(this.tissues==tiss)]))
  
  mtext(label.ages[names(pos.ages)], side=3, at=pos.ages, line=0, cex=0.55)
}


index=1:(length(this.samples)-1)
diff.species=index[unlist(lapply(index, function(x) this.species[x]!=this.species[x+1]))]
diff.stages=index[unlist(lapply(index, function(x) this.ages[x]!=this.ages[x+1]))]
diff.tissues=index[unlist(lapply(index, function(x) this.tissues[x]!=this.tissues[x+1]))]


## color rectangles for diff.species


width=diff(xpos)[1]

col.rect=rep(NA, length(this.samples))
col.rect[which(this.species=="Mouse")]="black"
col.rect[which(this.species=="Rat")]="gray40"
col.rect[which(this.species=="Chicken")]="gray80"
col.rect[which(this.samples=="Mouse_Testis_EarlyEmbryo")]=NA

rect(xpos-width/2, -0.019, xpos+width/2, -0.01, col=col.rect, border=NA, xpd=NA)

smally=abs(diff(ypos)[1])/2
## smally=0

for(i in setdiff(diff.stages, diff.tissues)){
  segments((xpos[i]+xpos[i+1])/2, 0-smally, (xpos[i]+xpos[i+1])/2, 1+smally, col="white", lwd=1)
}


for(i in diff.tissues){
  segments((xpos[i]+xpos[i+1])/2, 0-smally, (xpos[i]+xpos[i+1])/2, 1+smally)
}


boundaries=unlist(lapply(c("Brain", "Kidney", "Liver","Testis", "ImmuneCells"), function(x) which(markers$Organ==x)[1]))
## immuneboundary=which(markers$Organ=="ImmuneCells")[1]
yw=abs(diff(ypos))[1]/2

segments(-0.1, ypos[boundaries]+yw,1.1,ypos[boundaries]+yw)

## mtext for gene names

mtext(markers$GeneName, at=ypos, side=4, line=0, cex=0.5, font=3, las=2)

## cell type labels

cell.types=markers$CellType
tissues=markers$Organ
tissues[which(tissues%in%c("BrainMyeloid", "BrainNonMyeloid"))]="Brain"

tisscell=paste(tissues, cell.types, sep=",")
midy=unlist(lapply(unique(tisscell), function(x) mean(ypos[which(tisscell==x)])))

tcell=unlist(lapply(unique(tisscell), function(x) unlist(strsplit(x, split=","))[2]))
tcell=unlist(lapply(tcell, function(x) paste(unlist(strsplit(x, split="_")), collapse=" ")))
tissc=unlist(lapply(unique(tisscell), function(x) unlist(strsplit(x, split=","))[1]))

this.colors=col.tissues[tissc]
this.colors[which(tissc=="ImmuneCells")]="gray20"

mtext(tcell, side=4, at=midy, cex=0.55, las=2, line=3.95, col=this.colors, font=2)

tinyy=abs(diff(ypos)[1])/10

for(tc in unique(tisscell)){
  this.y=ypos[which(tisscell==tc)]
  if(length(this.y)>=1){
    t=unlist(strsplit(tc, split=","))[1]

    if(t%in%names(col.tissues)){
       segments(1.11, min(this.y)-tinyy, 1.11, max(this.y)+tinyy, col=col.tissues[t], xpd=NA)
     } else{
       segments(1.11, min(this.y)-tinyy, 1.11, max(this.y)+tinyy, col="gray20", xpd=NA)
     }
  }
}


mtext("species", side=4, at=-0.0135, cex=0.5, las=2, line=0.25)

#####################################################################################

## legend for the z-score

par(mar=c(2.1,1.5,1.5,45))

min.val=min(as.numeric(exp.zscore),na.rm=T)
max.val=max(as.numeric(exp.zscore), na.rm=T)
zlim=c(min.val, max.val)

z <- seq(min.val, max.val, length = 50)
## xax=pretty(z, ne=3)
xax=c(-2,0,2,4)
xax=xax[which(xax>=min.val & xax<=max.val)]

image(x=z, z = matrix(z, ncol = 1), col = topo.colors(50), zlim=zlim, xlim=range(xax), xaxt="n" ,yaxt="n")

par(tck=-0.25)
axis(side=1, at = xax, labels = xax,cex.axis=0.85, mgp=c(2,0.15,0))

mtext("expression z-score\n(log2-transformed TPM)",side=4, at=0.0, cex=0.65,line=0.5, las=2)


legend("topleft", legend=c("mouse", "rat"), fill=c("black", "gray40"), border=c("black", "gray40"), bty="n", inset=c(2.75,-0.5), xpd=NA, cex=0.95, horiz=T)

#####################################################################################

dev.off()

#######################################################################################

