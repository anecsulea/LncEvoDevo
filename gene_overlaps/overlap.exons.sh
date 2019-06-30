#!/bin/bash

export sp=$1

#####################################################################

export path=LncEvoDevo

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathResults=${path}/results/gene_overlaps/${sp}
export pathScripts=${path}/scripts/gene_overlaps

export release=94

#####################################################################

perl ${pathScripts}/overlap.exons.pl --pathExonBlocks1=${pathStringTie}/ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}.txt --pathExonBlocks2=${pathEnsembl}/ExonBlocks_AllTranscripts_Ensembl${release}.txt --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --biotypes="miRNA" --pathOutput=${pathResults}/ExonOverlap_StringTie_Ensembl_miRNA.txt 

#####################################################################
