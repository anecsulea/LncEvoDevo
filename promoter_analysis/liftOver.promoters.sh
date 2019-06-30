#!/bin/bash

export ref=$1
export tg=$2

##################################################################
##################################################################

export path=LncEvoDevo
export pathPromoters=${path}/results/promoter_analysis
export pathAlignments=${path}/data/genome_alignments

export release=94

##################################################################

if [ ${ref} = "Mouse" ]&&[ ${tg} = "Rat" ]; then
    export pathLiftOver=${pathAlignments}/mm10ToRn6.over.chain.gz
fi

##################################################################

if [ ${ref} = "Rat" ]&&[ ${tg} = "Mouse" ]; then
    export pathLiftOver=${pathAlignments}/rn6ToMm10.over.chain.gz
fi

##################################################################

liftOver -minMatch=0.5 ${pathPromoters}/${ref}/PromoterCoords_1kb_FilteredTranscripts_StringTie_Ensembl${release}.bed ${pathLiftOver} ${pathPromoters}/${tg}/From${ref}_PromoterCoords_1kb_FilteredTranscripts_StringTie_Ensembl${release}.bed ${pathPromoters}/${tg}/From${ref}_PromoterCoords_1kb_FilteredTranscripts_StringTie_Ensembl${release}.unmapped

##################################################################
