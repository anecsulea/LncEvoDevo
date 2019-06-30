##################################################################################

path="LncEvoDevo/"
pathResults=paste(path, "results/mutual_information_network/lncRNAs_only/", sep="")
pathTables=paste(path, "supplementary_tables/", sep="")

##################################################################################

go=read.table(paste(pathResults, "GOEnrichment_BiologicalProcess_lncRNAs_min50PCConnections.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

go=go[,c("ID.Mouse", "ID.Rat", "NbPCConnections", "GOAccession", "GOName", "Enrichment", "FDR")]

colnames(go)=c("ID.Mouse", "ID.Rat", "NbConnectedGenes",  "GOAccession", "GOName", "Enrichment", "FDR")

write.table(go, file=paste(pathTables, file="SupplementaryTable8.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

##################################################################################


