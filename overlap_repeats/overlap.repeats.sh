#!/bin/bash

export sp=$1
export annot=$2

##############################################################

export path=LncEvoDevo

##############################################################

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathRepeatMasker=${path}/data/RepeatMasker/${sp}
export pathResults=${path}/results/overlap_repeats/${sp}
export pathScripts=${path}/scripts/overlap_repeats

export release=94

##############################################################

if [ ${annot} = "Ensembl" ]; then
    export pathGTF=${pathEnsembl}/FilteredTranscripts_Ensembl${release}.gtf
    export prefixOutput=FilteredTranscripts_Ensembl${release}
fi

##############################################################

if [ ${annot} = "StringTie" ]; then
    export pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf
    export prefixOutput=FilteredTranscripts_StringTie_Ensembl${release}
fi

##############################################################

perl ${pathScripts}/overlap.repeats.pl --pathAnnotGTF=${pathGTF} --pathRepeatMasker=${pathRepeatMasker}/RepeatMasker_UCSC.txt.gz --pathOutputSense=${pathResults}/OverlapRepeats_SenseStrand_${prefixOutput}.txt --pathOutputAntisense=${pathResults}/OverlapRepeats_AntisenseStrand_${prefixOutput}.txt --pathOutputBothStrands=${pathResults}/OverlapRepeats_BothStrands_${prefixOutput}.txt

##############################################################
