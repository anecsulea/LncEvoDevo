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

perl ${pathScripts}/extract.random.intergenic.regions.pl --pathIntergenicRegions=${pathResults}/IntergenicRegions_${prefix}.txt --regionSize=5000 --minFlankingDistance=5000 --pathOutput=${pathResults}/ResampledIntergenicRegions_Size5kb_MinDistance5kb_${prefix}.txt

#############################################################################
