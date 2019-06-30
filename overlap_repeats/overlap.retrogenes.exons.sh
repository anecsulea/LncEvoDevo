#!/bin/bash

export sp=$1
export annot=$2

##############################################################

export path=LncEvoDevo
    
##############################################################

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathRetroCarelli=${path}/data/Retrogenes_Carelli2016/${sp}
export pathUCSC=${path}/data/UCSC_annotations/${sp}
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

if [ ${sp} = "Mouse" ]; then
    export suffixRetrogenesCarelli=Retrogenes_mm10.bed
    export suffixRetrogenesUCSC=RetroGenes_V6_mm10_Coords.bed
    export pathRetrogenes=${pathRetroCarelli}/${suffixRetrogenesCarelli},${pathUCSC}/${suffixRetrogenesUCSC}
fi

if [ ${sp} = "Rat" ]; then
    export suffixRetrogenesCarelli=Retrogenes_rn6.bed
    export pathRetrogenes=${pathRetroCarelli}/${suffixRetrogenesCarelli}
fi

if [ ${sp} = "Chicken" ]; then
    export suffixRetrogenesCarelli=Retrogenes_galGal5.bed
    export pathRetrogenes=${pathRetroCarelli}/${suffixRetrogenesCarelli}
fi

##############################################################

perl ${pathScripts}/overlap.retrogenes.exons.pl --pathExonBlocks=${pathExonBlocks} --pathRetrogenes=${pathRetrogenes} --type=both  --pathOutput=${pathResults}/OverlapRetrogenes_BothStrands_${prefixOutput}.txt


perl ${pathScripts}/overlap.retrogenes.exons.pl --pathExonBlocks=${pathExonBlocks} --pathRetrogenes=${pathRetrogenes} --type=sense --pathOutput=${pathResults}/OverlapRetrogenes_SenseStrand_${prefixOutput}.txt

perl ${pathScripts}/overlap.retrogenes.exons.pl --pathExonBlocks=${pathExonBlocks} --pathRetrogenes=${pathRetrogenes} --type=antisense --pathOutput=${pathResults}/OverlapRetrogenes_AntisenseStrand_${prefixOutput}.txt


##############################################################
