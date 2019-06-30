#!/bin/bash

export species=$1

#############################################################################

export path=LncEvoDevo
    
export pathStringTie=${path}/results/stringtie_assembly/${species}/combined
export pathResults=${path}/results/intergenic_regions/${species}
export pathScripts=${path}/scripts/intergenic_regions

export release=94
export prefix=FilteredTranscripts_StringTie_Ensembl${release}

#############################################################################

if [ -e ${pathResults} ]; then
    echo "path out already there"
else
    mkdir ${pathResults}
fi

#############################################################################

perl ${pathScripts}/extract.intergenic.regions.pl --pathExonBlocks=${pathStringTie}/ExonBlocks_${prefix}.txt --pathOutput=${pathResults}/IntergenicRegions_${prefix}.txt

#############################################################################

