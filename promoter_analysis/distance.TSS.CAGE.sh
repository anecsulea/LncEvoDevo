#!/bin/bash

export sp=$1

##########################################################################

export path=LncEvoDevo
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathData=${path}/data/CAGE
export pathResults=${path}/results/promoter_analysis/${sp}
export pathScripts=${path}/scripts/promoter_analysis

##########################################################################

export release=94
export prefix=FilteredTranscripts_StringTie_Ensembl${release}

if [ ${sp} = "Mouse" ]; then
    export pathCAGE=${pathData}/mm10.cage_peak_phase1and2combined_coord.bed.gz
fi

if [ ${sp} = "Rat" ]; then
    export pathCAGE=${pathData}/rn6.cage_peak_coord.bed.gz
fi

##########################################################################

perl ${pathScripts}/distance.TSS.CAGE.pl --pathGTF=${pathStringTie}/${prefix}.gtf --pathCAGE=${pathCAGE} --pathOutput=${pathResults}/DistanceCAGE_${prefix}.txt

##########################################################################
