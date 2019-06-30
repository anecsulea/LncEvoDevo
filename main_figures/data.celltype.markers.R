################################################################################

source("parameters.R")

pathTables=paste(path, "supplementary_tables/", sep="")

################################################################################

markers=read.table(paste(pathTables, "SupplementaryTable3.txt", sep=""), h=T, sep="\t", quote="\"", stringsAsFactors=F)

save(markers, file="RData/data.celltype.markers.RData")

################################################################################


