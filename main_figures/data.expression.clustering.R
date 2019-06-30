########################################################################

library(ape)


maxFDR=0.01

########################################################################

load("RData/data.expression.ortho.RData")

########################################################################

tpm.orthopc=exportho.mr[which(exportho.mr$GeneType=="protein_coding"), which(!colnames(exportho.mr)%in%c("ID", "GeneType"))]
tpm.ortholnc=exportho.mr[which(exportho.mr$GeneType=="lncRNA"), which(!colnames(exportho.mr)%in%c("ID", "GeneType"))]

samples=colnames(tpm.orthopc)
species=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))
tissues=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[2]))
ages=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[3]))
ages=unlist(lapply(ages, function(x) substr(x, 1, nchar(x)-1)))

tpm.orthopc=tpm.orthopc[,samples]
tpm.ortholnc=tpm.ortholnc[,samples]


idmouse.ortholnc=unlist(lapply(rownames(tpm.ortholnc), function(x) unlist(strsplit(x, split="_"))[1]))
idrat.ortholnc=unlist(lapply(rownames(tpm.ortholnc), function(x) unlist(strsplit(x, split="_"))[2]))

idmouse.orthopc=unlist(lapply(rownames(tpm.orthopc), function(x) unlist(strsplit(x, split="_"))[1]))
idrat.orthopc=unlist(lapply(rownames(tpm.orthopc), function(x) unlist(strsplit(x, split="_"))[2]))



########################################################################

cormat.pc.tpm=cor(log2(tpm.orthopc+1), method="spearman", use="complete.obs")

cormat.lnc.tpm=cor(log2(tpm.ortholnc+1), method="spearman", use="complete.obs")

########################################################################

dist.pc.tpm=as.dist(1-cormat.pc.tpm)

dist.lnc.tpm=as.dist(1-cormat.lnc.tpm)

hcl.pc.tpm=hclust(dist.pc.tpm)

hcl.lnc.tpm=hclust(dist.lnc.tpm)

sample.order.pc.tpm=order.dendrogram(as.dendrogram(hcl.pc.tpm))

sample.order.lnc.tpm=order.dendrogram(as.dendrogram(hcl.lnc.tpm))

tree.pc.tpm=as.phylo(hcl.pc.tpm)

tree.lnc.tpm=as.phylo(hcl.lnc.tpm)

########################################################################

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

#############################################################################

cormat.pc.tpm.diffexp=cor(log2(tpm.orthopc.diffexp+1), method="spearman", use="complete.obs")

cormat.lnc.tpm.diffexp=cor(log2(tpm.ortholnc.diffexp+1), method="spearman", use="complete.obs")

########################################################################

dist.pc.tpm.diffexp=as.dist(1-cormat.pc.tpm.diffexp)

dist.lnc.tpm.diffexp=as.dist(1-cormat.lnc.tpm.diffexp)

hcl.pc.tpm.diffexp=hclust(dist.pc.tpm.diffexp)

hcl.lnc.tpm.diffexp=hclust(dist.lnc.tpm.diffexp)

sample.order.pc.tpm.diffexp=order.dendrogram(as.dendrogram(hcl.pc.tpm.diffexp))

sample.order.lnc.tpm.diffexp=order.dendrogram(as.dendrogram(hcl.lnc.tpm.diffexp))

tree.pc.tpm.diffexp=as.phylo(hcl.pc.tpm.diffexp)

tree.lnc.tpm.diffexp=as.phylo(hcl.lnc.tpm.diffexp)

#############################################################################

save(list=c("cormat.pc.tpm", "cormat.lnc.tpm", "sample.order.pc.tpm", "sample.order.lnc.tpm", "tree.pc.tpm", "tree.lnc.tpm","cormat.pc.tpm.diffexp", "cormat.lnc.tpm.diffexp", "sample.order.pc.tpm.diffexp", "sample.order.lnc.tpm.diffexp", "tree.pc.tpm.diffexp", "tree.lnc.tpm.diffexp"), file="RData/data.expression.clustering.RData")

########################################################################
