########################################################################

library(ape)

load("RData/data.developmental.TF.RData")

########################################################################

load("RData/data.expression.ortho.RData")

########################################################################

tpm.orthopc=exportho.mrc[which(exportho.mrc$GeneType=="protein_coding"), which(!colnames(exportho.mrc)%in%c("ID", "GeneType"))]

cormat.pc.tpm=cor(log2(tpm.orthopc+1), method="spearman", use="complete.obs")

dist.pc.tpm=as.dist(1-cormat.pc.tpm)

hcl.pc.tpm=hclust(dist.pc.tpm)

sample.order.pc.tpm=order.dendrogram(as.dendrogram(hcl.pc.tpm))

tree.pc.tpm=as.phylo(hcl.pc.tpm)

########################################################################

id.mouse=unlist(lapply(rownames(tpm.orthopc), function(x) unlist(strsplit(x, split="_"))[1]))

tpm.orthopc.devTF=tpm.orthopc[which(id.mouse%in%dev.tf),]

cormat.devTF.tpm=cor(log2(tpm.orthopc.devTF+1), method="spearman", use="complete.obs")

dist.devTF.tpm=as.dist(1-cormat.devTF.tpm)

hcl.devTF.tpm=hclust(dist.devTF.tpm)

sample.order.devTF.tpm=order.dendrogram(as.dendrogram(hcl.devTF.tpm))

tree.devTF.tpm=as.phylo(hcl.devTF.tpm)

########################################################################

save(list=c("cormat.pc.tpm", "sample.order.pc.tpm", "tree.pc.tpm", "cormat.devTF.tpm", "sample.order.devTF.tpm", "tree.devTF.tpm"), file="RData/data.expression.clustering.mrc.RData")

########################################################################
