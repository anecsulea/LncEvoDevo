#!/bin/bash

export sp=$1
export annot=$2

######################################################################

export path=LncEvoDevo

######################################################################

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined

export release=94

######################################################################

if [ ${annot} = "Ensembl" ]; then
    if [ -e ${pathEnsembl}/bed_parts ]; then
	echo "output dir already there"
    else
	mkdir ${pathEnsembl}/bed_parts 
    fi

    split -l 10000 ${pathEnsembl}/ExonBlocks_FilteredTranscripts_Ensembl${release}.bed ${pathEnsembl}/bed_parts/ExonBlocks_FilteredTranscripts_Ensembl${release}_part  ## --additional-suffix=.bed
    
    for file in ${pathEnsembl}/bed_parts/*
    do
	mv "$file" "$file.bed"
    done
    
fi

######################################################################

if [ ${annot} = "StringTie" ]; then
   if [ -e ${pathStringTie}/bed_parts ]; then
	echo "output dir already there"
    else
	mkdir ${pathStringTie}/bed_parts 
    fi

    split -l 10000 ${pathStringTie}/ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}.bed ${pathStringTie}/bed_parts/ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}_part # --additional-suffix=.bed

    for file in ${pathStringTie}/bed_parts/*
    do
	mv "$file" "$file.bed"
    done
fi

######################################################################
