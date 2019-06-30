#!/bin/bash

export sp=$1

#########################################################################

export path=LncEvoDevo

export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathResults=${path}/results/gene_overlaps/${sp}
export pathScripts=${path}/scripts/gene_overlaps

export release=94

#########################################################################

perl ${pathScripts}/classify.gene.regions.pl --pathAllExonBlocks=${pathStringTie}/ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}.txt --pathSelectedExonBlocks=${pathStringTie}/ExonBlocks_ProteinCodingGenes_FilteredTranscripts_StringTie_Ensembl${release}.txt --pathSelectedIntronBlocks=${pathStringTie}/IntronBlocks_ProteinCodingGenes_FilteredTranscripts_StringTie_Ensembl${release}.txt --pathOutput=${pathResults}/GeneClassification_ProteinCodingGenes_FilteredTranscripts_StringTie_Ensembl${release}.txt


perl ${pathScripts}/classify.gene.regions.pl --pathAllExonBlocks=${pathStringTie}/ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}.txt --pathSelectedExonBlocks=${pathStringTie}/ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}.txt --pathSelectedIntronBlocks=${pathStringTie}/IntronBlocks_AllGenes_FilteredTranscripts_StringTie_Ensembl${release}.txt --pathOutput=${pathResults}/GeneClassification_AllGenes_FilteredTranscripts_StringTie_Ensembl${release}.txt


#########################################################################







