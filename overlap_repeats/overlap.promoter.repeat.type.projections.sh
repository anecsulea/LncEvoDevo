#!/bin/bash

export ref=$1
export tg=$2

##############################################################
    
export path=LncEvoDevo

##############################################################

export pathRepeatMasker=${path}/data/RepeatMasker/${tg}
export pathProjections=${path}/results/exon_projections
export pathResults=${path}/results/overlap_repeats/${tg}
export pathScripts=${path}/scripts/overlap_repeats

export release=94

##############################################################

export pathGTF=${pathProjections}/From${ref}_To${tg}_ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}_FilteredProjectedExons_Step2.gtf
export prefixOutput=From${ref}_To${tg}_ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}_FilteredProjectedExons_Step2

##############################################################

perl ${pathScripts}/overlap.promoter.repeat.type.pl --pathGTF=${pathGTF} --promoterSize=1000 --pathRepeatMasker=${pathRepeatMasker}/RepeatMasker_UCSC.txt.gz --type=repFamily --repeatClasses=any --pathOutput=${pathResults}/OverlapRepeats_TEFamily_BothStrands_Promoters1kb_${prefixOutput}.txt

perl ${pathScripts}/overlap.promoter.repeat.type.pl --pathGTF=${pathGTF} --promoterSize=1000 --pathRepeatMasker=${pathRepeatMasker}/RepeatMasker_UCSC.txt.gz --type=repClass --repeatClasses=any --pathOutput=${pathResults}/OverlapRepeats_TEClass_BothStrands_Promoters1kb_${prefixOutput}.txt

##############################################################

perl ${pathScripts}/overlap.promoter.repeat.type.pl --pathGTF=${pathGTF} --promoterSize=5000 --pathRepeatMasker=${pathRepeatMasker}/RepeatMasker_UCSC.txt.gz --type=repFamily --repeatClasses=any --pathOutput=${pathResults}/OverlapRepeats_TEFamily_BothStrands_Promoters5kb_${prefixOutput}.txt

perl ${pathScripts}/overlap.promoter.repeat.type.pl --pathGTF=${pathGTF} --promoterSize=5000 --pathRepeatMasker=${pathRepeatMasker}/RepeatMasker_UCSC.txt.gz --type=repClass --repeatClasses=any --pathOutput=${pathResults}/OverlapRepeats_TEClass_BothStrands_Promoters5kb_${prefixOutput}.txt

##############################################################
