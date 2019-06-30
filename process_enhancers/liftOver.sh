#!/bin/bash

#########################################################

export path=LncEvoDevo
export pathEnhancers=${path}/data/enhancers/Encode_YueLab/predicted_enhancer_mouse
export pathLiftOver=${path}/data/genome_alignments

#########################################################

liftOver ${pathEnhancers}/all_enhancers_500bp_mm9.bed ${pathLiftOver}/mm9ToMm10.over.chain.gz ${pathEnhancers}/all_enhancers_500bp_mm10.bed ${pathEnhancers}/all_enhancers_500bp_mm9.unmapped 

#########################################################

