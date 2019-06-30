#############################################################################

load("RData/data.expression.ortho.RData")

library(ade4)

maxFDR=0.01

#############################################################################

tpm.orthopc=exportho.mr[which(exportho.mr$GeneType=="protein_coding"), which(!colnames(exportho.mr)%in%c("ID", "GeneType"))]
tpm.ortholnc=exportho.mr[which(exportho.mr$GeneType=="lncRNA"), which(!colnames(exportho.mr)%in%c("ID", "GeneType"))]

hasnapc=apply(tpm.orthopc, 1, function(x) any(is.na(x)))
tpm.orthopc=tpm.orthopc[which(!hasnapc),]

hasnalnc=apply(tpm.ortholnc, 1, function(x) any(is.na(x)))
tpm.ortholnc=tpm.ortholnc[which(!hasnalnc),]

idmouse.ortholnc=unlist(lapply(rownames(tpm.ortholnc), function(x) unlist(strsplit(x, split="_"))[1]))
idrat.ortholnc=unlist(lapply(rownames(tpm.ortholnc), function(x) unlist(strsplit(x, split="_"))[2]))

idmouse.orthopc=unlist(lapply(rownames(tpm.orthopc), function(x) unlist(strsplit(x, split="_"))[1]))
idrat.orthopc=unlist(lapply(rownames(tpm.orthopc), function(x) unlist(strsplit(x, split="_"))[2]))


#############################################################################

samples=colnames(tpm.orthopc)
species=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))
tissues=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[2]))
ages=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[3]))
ages=unlist(lapply(ages, function(x) substr(x, 1, nchar(x)-1)))

tpm.ortholnc=tpm.ortholnc[,samples]


pca.pc=dudi.pca(log2(tpm.orthopc+1), center=T, scale=T, scannf=F, nf=5)
pca.lnc=dudi.pca(log2(tpm.ortholnc+1), center=T, scale=T, scannf=F, nf=5)

#############################################################################

## extract differentially expressed protein-coding and lncRNAs

load("RData/data.diffexp.Mouse.RData")
de.global.mouse=de.global.allreads

diffexp.mouse=de.global.mouse$GeneID[which(apply(de.global.mouse[,grep("FDR", colnames(de.global.mouse))], 1, min)<maxFDR)]


load("RData/data.diffexp.Rat.RData")
de.global.rat=de.global.allreads

diffexp.rat=de.global.rat$GeneID[which(apply(de.global.rat[,grep("FDR", colnames(de.global.rat))], 1, min)<maxFDR)]

#############################################################################

tpm.ortholnc.diffexp=tpm.ortholnc[which(idmouse.ortholnc%in%diffexp.mouse | idrat.ortholnc%in%diffexp.rat),]
tpm.orthopc.diffexp=tpm.orthopc[which(idmouse.orthopc%in%diffexp.mouse | idrat.orthopc%in%diffexp.rat),]

pca.lnc.diffexp=dudi.pca(log2(tpm.ortholnc.diffexp+1), center=T, scale=T, scannf=F, nf=5)
pca.pc.diffexp=dudi.pca(log2(tpm.orthopc.diffexp+1), center=T, scale=T, scannf=F, nf=5)

#############################################################################

save(list=c("pca.pc", "pca.lnc", "pca.lnc.diffexp", "pca.pc.diffexp"), file="RData/data.pca.RData")

#############################################################################
