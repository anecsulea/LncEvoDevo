#!/bin/bash

export sp=$1

############################################################################################

export path=LncEvoDevo
export pathAnnot=${path}/data/ensembl_annotations/${sp}
export pathResults=${path}/results/CSF/intron_regions/${sp}
export pathScripts=${path}/scripts/coding_potential/CSF

export release=84

############################################################################################

if [ -e ${pathResults} ]; then
    echo "path results exists"
else
    mkdir ${pathResults}
fi

perl ${pathScripts}/extract.introns.pl --pathExonBlocks=${pathAnnot}/ExonBlocks_Ensembl${release}_AllTranscripts.txt --minIntronSize=500 --maxIntronSize=10000 --excludeLength=200 --pathOutput=${pathResults}/IntronRegions_MinSize500_MaxSize10000_ExcludeLength200_Ensembl${release}.txt 

############################################################################################
 
