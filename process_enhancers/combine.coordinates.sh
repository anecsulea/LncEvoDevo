#!/bin/bash

##############################################################

export path=LncEvoDevo
export pathData=${path}/data/enhancers/Encode_YueLab/predicted_enhancer_mouse
export pathScripts=${path}/scripts/process_enhancers

##############################################################

export paths=""

for file in `ls ${pathData} | grep final`
do
    export paths=${pathData}/${file},${paths}
done

echo ${paths}

##############################################################

perl ${pathScripts}/combine.coordinates.pl --pathsEnhancers=${paths} --halfSize=250 --pathOutput=${pathData}/all_enhancers_500bp_mm9.bed

##############################################################
