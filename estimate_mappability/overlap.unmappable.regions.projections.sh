#!/bin/bash

export sp1=$1
export sp2=$2
export annot=$3

####################################################################################

export path=LncEvoDevo

####################################################################################

export pathStringTie=${path}/results/stringtie_assembly
export pathEnsembl=${path}/data/ensembl_annotations
export pathUCSC=${path}/data/UCSC_sequences
export pathMappability=${path}/results/mappability/${sp2}
export pathProjections=${path}/results/exon_projections
export pathResults=${path}/results/overlap_unmappable_regions/${sp1}
export pathScripts=${path}/scripts/estimate_mappability

export release=89

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

if [ ${sp1} = ${sp2} ]; then
    echo "not looking at ${sp1} against itself"
    exit
fi 

####################################################################################

perl ${pathScripts}/overlap.unmappable.regions.projections.pl --pathProjectedExons=${pathProjections}/From${sp1}_To${sp2}_${prefix}_FilteredProjectedExons_Step2.txt --pathUnmappableRegions=${pathMappability}/unmappable_regions.txt --pathOutput=${pathResults}/OverlapUnmappableRegions_From${sp1}_To${sp2}_${prefix}_FilteredProjectedExons_Step2.txt

####################################################################################
