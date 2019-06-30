############################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathResults=paste(path, "results/mutual_information_network/", sep="")

options(scipen=999)

############################################################################

## read expression data

tpm=read.table(paste(pathDatasets, "SupplementaryDataset6/NormalizedTPM_OrthoGenes_Mouse_Rat.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

rownames(tpm)=tpm$ID

info=tpm[,c(1,2)]
tpm=tpm[,-c(1,2)]

############################################################################

## maximum expression 1
## expression above 0 in at least 5 samples

minexp=1
minsamples=5

maxexp.mouse=apply(tpm[,grep("Mouse", colnames(tpm))], 1, max)
maxexp.rat=apply(tpm[,grep("Rat", colnames(tpm))], 1, max)

nbsamples.mouse=apply(tpm[,grep("Mouse", colnames(tpm))], 1, function(x) length(which(x>0)))
nbsamples.rat=apply(tpm[,grep("Rat", colnames(tpm))], 1, function(x) length(which(x>0)))

selected=which(nbsamples.mouse>=minsamples & nbsamples.rat>=minsamples & maxexp.mouse>=minexp & maxexp.rat>=minexp)

tpm=tpm[selected,]
info=info[selected,]

############################################################################

tpm.mouse=tpm[,grep("Mouse", colnames(tpm))]
tpm.mouse$Gene=info$ID
tpm.mouse=tpm.mouse[,c("Gene", setdiff(colnames(tpm.mouse), "Gene"))]

tpm.rat=tpm[,grep("Rat", colnames(tpm))]
tpm.rat$Gene=info$ID
tpm.rat=tpm.rat[,c("Gene", setdiff(colnames(tpm.rat), "Gene"))]

############################################################################

write.table(tpm.mouse, file=paste(pathResults, "lncRNAs_only/Mouse/TPM.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
write.table(tpm.rat, file=paste(pathResults, "lncRNAs_only/Rat/TPM.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
writeLines(info$ID[which(info$GeneType=="lncRNA")], con=paste(pathResults, "lncRNAs_only/LncRNAs.txt", sep=""))

write.table(tpm.mouse, file=paste(pathResults, "all_genes/Mouse/TPM.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
write.table(tpm.rat, file=paste(pathResults, "all_genes/Rat/TPM.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
writeLines(info$ID, con=paste(pathResults, "all_genes/AllGenes.txt", sep=""))

############################################################################
