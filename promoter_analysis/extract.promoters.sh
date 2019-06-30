#!/bin/bash

export sp=$1
export annot=$2

##############################################################################

export path=LncEvoDevo
export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathFasta=${path}/data/genome_indexes/${sp}
export pathUCSC=${path}/data/UCSC_sequences/${sp}
export pathResults=${path}/results/promoter_analysis/${sp}
export pathScripts=${path}/scripts/promoter_analysis

##############################################################################

export release=94
export genomerelease=87

##############################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathExons=${pathEnsembl}
    export suffixExons=FilteredTranscripts_Ensembl${release}
fi

if [ ${annot} = "StringTie" ]; then
    export pathExons=${pathStringTie}
    export suffixExons=FilteredTranscripts_StringTie_Ensembl${release}
fi

#############################################################################

perl ${pathScripts}/extract.promoters.pl --pathGTF=${pathExons}/${suffixExons}.gtf  --pathGenomeSequence=${pathFasta}/genome_ensembl${genomerelease}.fa --promoterSize=1000 --pathChromosomeCorrespondence=${pathUCSC}/chromosomes_Ensembl_UCSC.txt --pathOutputCoords=${pathResults}/PromoterCoords_1kb_${suffixExons}.bed  --pathOutputSequences=${pathResults}/PromoterSequences_1kb_${suffixExons}.fa

#############################################################################

 perl ${pathScripts}/extract.promoters.pl --pathGTF=${pathExons}/${suffixExons}.gtf  --pathGenomeSequence=${pathFasta}/genome_ensembl${genomerelease}.fa --promoterSize=400 --pathChromosomeCorrespondence=${pathUCSC}/chromosomes_Ensembl_UCSC.txt --pathOutputCoords=${pathResults}/PromoterCoords_400bp_${suffixExons}.bed  --pathOutputSequences=${pathResults}/PromoterSequences_400bp_${suffixExons}.fa

#############################################################################
