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

## overlap sense miRNA

perl ${pathScripts}/overlap.genes.pl --pathGTF1=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathGTF2=${pathEnsembl}/AllTranscripts_Ensembl${release}.gtf --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --biotypes="miRNA" --windowSizes=0,5000,10000,20000,50000,100000  --forbiddenTranscripts=ENSMUST00000127091,ENSMUST00000127664 --overlapSense="sense" --pathOutput=${pathResults}/GeneOverlap_Sense_StringTie_Ensembl_miRNA.txt 

#####################################################################

## overlap sense tRNA

perl ${pathScripts}/overlap.genes.pl --pathGTF1=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathGTF2=${pathEnsembl}/AllTranscripts_Ensembl${release}.gtf --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --biotypes="tRNA" --windowSizes=0,5000,10000,20000,50000,100000  --forbiddenTranscripts=ENSMUST00000127091,ENSMUST00000127664 --overlapSense="sense" --pathOutput=${pathResults}/GeneOverlap_Sense_StringTie_Ensembl_tRNA.txt 

#####################################################################
