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


perl ${pathScripts}/overlap.promoter.repeat.type.pl --pathGTF=${pathGTF} --promoterSize=1000 --pathRepeatMasker=${pathRepeatMasker}/RepeatMasker_UCSC.txt.gz --type=repFamily --repeatClasses=any --pathOutput=${pathResults}/OverlapRepeats_TEFamily_BothStrands_Promoters1kb_${prefixOutput}.txt

perl ${pathScripts}/overlap.promoter.repeat.type.pl --pathGTF=${pathGTF} --promoterSize=1000 --pathRepeatMasker=${pathRepeatMasker}/RepeatMasker_UCSC.txt.gz --type=repClass --repeatClasses=any --pathOutput=${pathResults}/OverlapRepeats_TEClass_BothStrands_Promoters1kb_${prefixOutput}.txt

##############################################################

perl ${pathScripts}/overlap.promoter.repeat.type.pl --pathGTF=${pathGTF} --promoterSize=5000 --pathRepeatMasker=${pathRepeatMasker}/RepeatMasker_UCSC.txt.gz --type=repFamily --repeatClasses=any --pathOutput=${pathResults}/OverlapRepeats_TEFamily_BothStrands_Promoters5kb_${prefixOutput}.txt

perl ${pathScripts}/overlap.promoter.repeat.type.pl --pathGTF=${pathGTF} --promoterSize=5000 --pathRepeatMasker=${pathRepeatMasker}/RepeatMasker_UCSC.txt.gz --type=repClass --repeatClasses=any --pathOutput=${pathResults}/OverlapRepeats_TEClass_BothStrands_Promoters5kb_${prefixOutput}.txt

##############################################################
