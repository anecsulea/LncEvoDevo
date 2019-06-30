#!/bin/bash

export sp=$1
export annot=$2

#####################################################################

export path=LncEvoDevo

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathResults=${path}/results/gene_overlaps/${sp}
export pathScripts=${path}/scripts/gene_overlaps

export release=94

#####################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathAnnot=${pathEnsembl}
    export filename=FilteredTranscripts_Ensembl${release}
fi 


if [ ${annot} = "StringTie" ]; then
    export pathAnnot=${pathStringTie}
    export filename=FilteredTranscripts_StringTie_Ensembl${release}
fi 

#####################################################################

perl ${pathScripts}/extract.intergenic.regions.pl --pathExonBlocks=${pathAnnot}/ExonBlocks_${filename}.txt --minDistance=5000 --minSize=1000 --pathOutput=${pathAnnot}/IntergenicRegions_MinDistance5kb_MinSize1kb_${filename}.txt 

#####################################################################
