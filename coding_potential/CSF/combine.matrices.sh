#!/bin/bash

export sp=$1
export aln=$2

############################################################################################

export path=LncEvoDevo
export pathMatrices=${path}/results/CSF/matrices/${aln}/${sp}
export pathScripts=${path}/scripts/coding_potential/CSF

############################################################################################

export pathsIntrons=""
export pathsCDS=""

for chr in {1..22} X Y Z W
do
    if [ -e ${pathMatrices}/matrices_introns_chr${chr}.txt ]; then
	export pathsIntrons=${pathMatrices}/matrices_introns_chr${chr}.txt,${pathsIntrons}
    fi

    if [ -e ${pathMatrices}/matrices_CDS_chr${chr}.txt ]; then
	export pathsCDS=${pathMatrices}/matrices_CDS_chr${chr}.txt,${pathsCDS}
    fi
done 

############################################################################################

perl ${pathScripts}/combine.matrices.pl --pathGeneticCode=${pathScripts}/standard_genetic_code.txt --pathsMatrices=${pathsIntrons} --pathOutput=${pathMatrices}/matrices_introns.txt

perl ${pathScripts}/combine.matrices.pl --pathGeneticCode=${pathScripts}/standard_genetic_code.txt --pathsMatrices=${pathsCDS} --pathOutput=${pathMatrices}/matrices_CDS.txt
    
############################################################################################
