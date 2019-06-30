#############################################################################

load("RData/data.expression.ortho.RData")

load("RData/data.developmental.TF.RData")

library(ade4)

#############################################################################

tpm.orthopc=exportho.mrc[which(exportho.mrc$GeneType=="protein_coding"), which(!colnames(exportho.mrc)%in%c("ID", "GeneType"))]

hasnapc=apply(tpm.orthopc, 1, function(x) any(is.na(x)))
tpm.orthopc=tpm.orthopc[which(!hasnapc),]

#############################################################################

pca.pc=dudi.pca(log2(tpm.orthopc+1), center=T, scale=T, scannf=F, nf=5)

#############################################################################

id.mouse=unlist(lapply(rownames(tpm.orthopc), function(x) unlist(strsplit(x, split="_"))[1]))
tpm.orthopc.devTF=tpm.orthopc[which(id.mouse%in%dev.tf),]

pca.pc.devTF=dudi.pca(log2(tpm.orthopc.devTF+1), center=T, scale=T, scannf=F, nf=5)

#############################################################################

save(list=c("pca.pc", "pca.pc.devTF"), file="RData/data.pca.mrc.RData")

#############################################################################
