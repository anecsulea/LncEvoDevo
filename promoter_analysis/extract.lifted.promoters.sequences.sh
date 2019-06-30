#!/bin/bash

export ref=$1
export tg=$2

##################################################################
##################################################################

export path=LncEvoDevo
export pathPromoters=${path}/results/promoter_analysis
export pathAlignments=${path}/data/genome_alignments
export pathFasta=${path}/data/genome_indexes
export pathScripts=${path}/scripts/promoter_analysis

export release=94
export genomerelease=87

##################################################################

perl ${pathScripts}/extract.lifted.promoters.sequences.pl --pathCoordinates=${pathPromoters}/${tg}/From${ref}_PromoterCoords_1kb_FilteredTranscripts_StringTie_Ensembl${release}.bed --outprefix=${tg} --pathGenomeSequence=${pathFasta}/${tg}/genome_ensembl${genomerelease}.fa --pathOutputSequences=${pathPromoters}/${tg}/From${ref}_PromoterCoords_1kb_FilteredTranscripts_StringTie_Ensembl${release}.fa

##################################################################
