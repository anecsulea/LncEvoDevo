#!/bin/bash

export species=$1
export type=$2

####################################################################################

export path=LncEvoDevo

####################################################################################

export pathHisat=${path}/results/hisat/${species}
export pathDocs=${path}/docs
export pathResults=${path}/results/kallisto_indexes/${species}
export pathAnnot=${path}/data/ensembl_annotations/${species}
export pathStringTie=${path}/results/stringtie_assembly/${species}/combined
export pathScripts=${path}/scripts/expression_estimation

export release=94

####################################################################################


if [ ${type} = "Ensembl" ]; then
    export pathFasta=${pathAnnot}/FilteredTranscripts_Ensembl${release}_noMT.fa
fi


if [ ${type} = "StringTie" ]; then
    export pathFasta=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}_noMT.fa
fi

####################################################################################

if [ -e ${pathResults} ]; then
    echo "path results already there"
else
    mkdir ${pathResults}
fi

####################################################################################

kallisto index -i ${pathResults}/${type} ${pathFasta}

####################################################################################
