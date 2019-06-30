#!/bin/bash

export sp=$1
export annot=$2


####################################################################################

export path=LncEvoDevo

####################################################################################

export pathStringTie=${path}/results/stringtie_assembly
export pathEnsembl=${path}/data/ensembl_annotations
export pathUCSC=${path}/data/UCSC_sequences
export pathMappability=${path}/results/mappability/${sp}
export pathResults=${path}/results/overlap_unmappable_regions/${sp}
export pathScripts=${path}/scripts/estimate_mappability

export release=94

####################################################################################

if [ ${annot} = "Ensembl" ]; then
    export prefix=FilteredTranscripts_Ensembl${release}
    export pathExons=${pathEnsembl}/${sp}
fi

if [ ${annot} = "StringTie" ]; then
    export prefix=FilteredTranscripts_StringTie_Ensembl${release}
    export pathExons=${pathStringTie}/${sp}/combined
fi

####################################################################################

perl ${pathScripts}/overlap.unmappable.regions.pl --pathExonBlocks=${pathExons}/ExonBlocks_${prefix}.txt --pathUnmappableRegions=${pathMappability}/unmappable_regions.txt --pathOutput=${pathResults}/OverlapUnmappableRegions_${prefix}.txt

####################################################################################
