#!/bin/bash

export sp="Mouse"

#########################################################################

export path=LncEvoDevo

export pathEnhancers=${path}/data/enhancers/Encode_YueLab/predicted_enhancer_mouse
export pathVista=${path}/data/enhancers/Vista
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathResults=${path}/results/overlap_enhancers/${sp}
export pathScripts=${path}/scripts/process_enhancers

export release=94

#########################################################################

perl ${pathScripts}/overlap.tss.enhancers.pl --pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathEnhancers=${pathEnhancers}/all_enhancers_500bp_mm10.bed --maxDistance=500 --pathOutput=${pathResults}/OverlapTSS_EncodeYueLabEnhancers_MaxDist500bp_FilteredTranscripts_StringTie_Ensembl${release}.txt


perl ${pathScripts}/overlap.tss.enhancers.pl --pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathEnhancers=${pathEnhancers}/all_enhancers_500bp_mm10.bed --maxDistance=0 --pathOutput=${pathResults}/OverlapTSS_EncodeYueLabEnhancers_MaxDist0bp_FilteredTranscripts_StringTie_Ensembl${release}.txt

perl ${pathScripts}/overlap.tss.enhancers.pl --pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathEnhancers=${pathEnhancers}/all_enhancers_500bp_mm10.bed --maxDistance=1000 --pathOutput=${pathResults}/OverlapTSS_EncodeYueLabEnhancers_MaxDist1kb_FilteredTranscripts_StringTie_Ensembl${release}.txt

#########################################################################

perl ${pathScripts}/overlap.tss.enhancers.pl --pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathEnhancers=${pathVista}/PositiveRegions_mm10.bed --maxDistance=0 --pathOutput=${pathResults}/OverlapTSS_Vista_MaxDist0bp_FilteredTranscripts_StringTie_Ensembl${release}.txt

perl ${pathScripts}/overlap.tss.enhancers.pl --pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathEnhancers=${pathVista}/PositiveRegions_mm10.bed --maxDistance=500 --pathOutput=${pathResults}/OverlapTSS_Vista_MaxDist500bp_FilteredTranscripts_StringTie_Ensembl${release}.txt


perl ${pathScripts}/overlap.tss.enhancers.pl --pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathEnhancers=${pathVista}/PositiveRegions_mm10.bed --maxDistance=1000 --pathOutput=${pathResults}/OverlapTSS_Vista_MaxDist1kb_FilteredTranscripts_StringTie_Ensembl${release}.txt

#########################################################################
