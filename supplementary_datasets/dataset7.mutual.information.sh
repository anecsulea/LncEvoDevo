#!/bin/bash

##########################################################################

export path=LncEvoDevo
export pathNetwork=${path}/results/mutual_information_network/lncRNAs_only
export pathResults=${path}/supplementary_datasets/SupplementaryDataset7

##########################################################################

for sp in Mouse Rat
do
    cp ${pathNetwork}/${sp}/network.txt ${pathResults}/FullNetwork_${sp}.txt
    gzip ${pathResults}/FullNetwork_${sp}.txt
    
    cp ${pathNetwork}/${sp}/filtered_network_FDR0.001.txt ${pathResults}/FilteredNetwork_FDR0.001_${sp}.txt
    gzip ${pathResults}/FilteredNetwork_FDR0.001_${sp}.txt

done

##########################################################################

cp ${pathNetwork}/common_interactions_FDR0.001.txt ${pathResults}/SharedInteractions_FDR0.001.txt
gzip ${pathResults}/SharedInteractions_FDR0.001.txt

cp ${pathNetwork}/GOEnrichment_BiologicalProcess_lncRNAs_min50PCConnections.txt ${pathResults}/GOEnrichment_BiologicalProcess_lncRNAs_Min50ConnectedGenes.txt
gzip ${pathResults}/GOEnrichment_BiologicalProcess_lncRNAs_Min50ConnectedGenes.txt

##########################################################################
