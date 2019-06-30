################################################################################

path="LncEvoDevo/"
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathResults=paste(path, "results/mutual_information_network/", sep="")
pathFigures=paste(path, "scripts/main_figures/",sep="")

################################################################################


options(scipen=999)

maxFDR=0.001 ## to define mutual information network

type="lncRNAs_only"

minConnections=50 ## minimum number of connected protein-coding genes to compute GO enrichment

maxGOFDR=0.1

############################################################################

GOresults=read.table(paste(pathResults, type, "/GOEnrichment_BiologicalProcess_lncRNAs_min", minConnections,"PCConnections.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote=F)
