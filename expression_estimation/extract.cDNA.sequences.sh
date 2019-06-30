#!/bin/bash

export sp=$1
export annot=$2

##############################################################

export path=LncEvoDevo

###########################################################################

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathBlastDB=${path}/data/blastn_databases/${sp}
export pathGenome=${path}/data/genome_indexes/${sp}
export pathScripts=${path}/scripts/expression_estimation

export release=94
export genomerelease=87

##############################################################

if [ ${annot} = "Ensembl" ]; then
    export pathGTF=${pathEnsembl}
    export suffix=FilteredTranscripts_Ensembl${release}
fi

if [ ${annot} = "StringTie" ]; then
    export pathGTF=${pathStringTie}
    export suffix=FilteredTranscripts_StringTie_Ensembl${release}
fi

##############################################################

perl ${pathScripts}/extract.cDNA.sequences.pl --pathAnnotGTF=${pathGTF}/${suffix}.gtf --forbiddenChromo=MT --pathGenomeSequence=${pathGenome}/genome_ensembl${genomerelease}.fa --pathOutput=${pathGTF}/${suffix}_noMT.fa

###############################################################
