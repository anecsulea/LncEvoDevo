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
    export pathExonBlocks=${pathEnsembl}/ExonBlocks_FilteredTranscripts_Ensembl${release}.txt
    export prefixOutput=FilteredTranscripts_Ensembl${release}
fi

##############################################################

if [ ${annot} = "StringTie" ]; then
    export pathExonBlocks=${pathStringTie}/ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}.txt
    export prefixOutput=FilteredTranscripts_StringTie_Ensembl${release}
fi

##############################################################

perl ${pathScripts}/overlap.exons.repeats.type.pl --pathExonBlocks=${pathExonBlocks} --pathRepeatMasker=${pathRepeatMasker}/RepeatMasker_UCSC.txt.gz --type=repFamily --repeatClasses=RNA,rRNA,tRNA,scRNA,snRNA,srpRNA --pathOutput=${pathResults}/OverlapRepeats_RNA_BothStrands_${prefixOutput}.txt

##############################################################

perl ${pathScripts}/overlap.exons.repeats.type.pl --pathExonBlocks=${pathExonBlocks} --pathRepeatMasker=${pathRepeatMasker}/RepeatMasker_UCSC.txt.gz --type=repFamily --repeatClasses=any --pathOutput=${pathResults}/OverlapRepeats_TEFamily_BothStrands_${prefixOutput}.txt

perl ${pathScripts}/overlap.exons.repeats.type.pl --pathExonBlocks=${pathExonBlocks} --pathRepeatMasker=${pathRepeatMasker}/RepeatMasker_UCSC.txt.gz --type=repClass --repeatClasses=any --pathOutput=${pathResults}/OverlapRepeats_TEClass_BothStrands_${prefixOutput}.txt

##############################################################