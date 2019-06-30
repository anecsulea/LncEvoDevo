#!/bin/bash


#########################################################

export path=LncEvoDevo

export pathEnhancers=${path}/data/enhancers/Encode_YueLab/predicted_enhancer_mouse
export pathAnnot=${path}/results/stringtie_assembly/Mouse/combined
export pathLiftOver=${path}/data/genome_alignments
export pathResults=${path}/results/overlap_enhancers/Mouse
export pathScripts=${path}/scripts/process_enhancers

export release=94

#########################################################

perl ${pathScripts}/overlap.exons.enhancers.pl --pathExonBlocks=${pathAnnot}/ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}.txt --pathEnhancers=${pathEnhancers}/all_enhancers_500bp_mm10.bed --pathOutput=${pathResults}/ExonicOverlap_EncodeYueLabEnhancers_ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}.txt 

#########################################################
