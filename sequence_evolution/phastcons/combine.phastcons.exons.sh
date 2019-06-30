#!/bin/bash

export species=$1
export annot=$2
export phast=$3

#####################################################################

export path=LncEvoDevo
    
export pathEnsembl=${path}/data/ensembl_annotations/${species}
export pathStringTie=${path}/results/stringtie_assembly/${species}/combined
export pathResults=${path}/results/sequence_evolution/phastcons/${species}/${phast}
export pathPhastCons=${path}/data/phastcons/${species}
export pathScripts=${path}/scripts/sequence_evolution/phastcons

export release=94

#####################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathExons=${pathEnsembl}
    export suffixExons=ExonBlocks_FilteredTranscripts_Ensembl${release}
fi

if [ ${annot} = "StringTie" ]; then
    export pathExons=${pathStringTie}
    export suffixExons=ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}
fi

#####################################################################

cp ${pathResults}/PhastCons_ExonAverage_chr1_${suffixExons}.txt ${pathResults}/PhastCons_ExonAverage_${suffixExons}.txt 
cp ${pathResults}/PhastCons_GeneAverage_chr1_${suffixExons}.txt ${pathResults}/PhastCons_GeneAverage_${suffixExons}.txt 

for chr in {2..19} X Y
do

    sed '1d'  ${pathResults}/PhastCons_ExonAverage_chr${chr}_${suffixExons}.txt >> ${pathResults}/PhastCons_ExonAverage_${suffixExons}.txt 

    sed '1d'  ${pathResults}/PhastCons_GeneAverage_chr${chr}_${suffixExons}.txt >> ${pathResults}/PhastCons_GeneAverage_${suffixExons}.txt 

done

#####################################################################
