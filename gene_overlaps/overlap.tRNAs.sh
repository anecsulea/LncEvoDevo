#!/bin/bash

export sp=$1
export annot=$2

##############################################################

export path=LncEvoDevo

##############################################################

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathUCSC=${path}/data/UCSC_annotations/${sp}
export pathResults=${path}/results/gene_overlaps/${sp}
export pathScripts=${path}/scripts/gene_overlaps

export release=94

##############################################################

if [ ${annot} = "Ensembl" ]; then
    export pathExonBlocks=${pathEnsembl}/ExonBlocks_FilteredTranscripts_Ensembl${release}.txt
    export prefixOutput=FilteredTranscripts_Ensembl${release}
fi

##############################################################

if [ ${annot} = "StringTie" ]; then
    export pathExonBlocks=${pathStringTie}/ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}.txt
    export prefixOutput=FilteredTranscripts_StringTie_Ensembl${release}
fi

##############################################################

perl ${pathScripts}/overlap.tRNAs.pl --pathExonBlocks=${pathExonBlocks} --pathTRNAs=${pathUCSC}/tRNAGenes.txt.gz --pathOutput=${pathResults}/Overlap_tRNAs_${prefixOutput}.txt

##############################################################
