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
export pathRepeatMasker=${path}/data/RepeatMasker/${sp}
export pathResults=${path}/results/gene_overlaps/${sp}
export pathScripts=${path}/scripts/gene_overlaps

export release=94

##############################################################

if [ ${annot} = "Ensembl" ]; then
    export pathAnnot=${pathEnsembl}
    export prefixAnnot=FilteredTranscripts_Ensembl${release}
fi

##############################################################

if [ ${annot} = "StringTie" ]; then
    export pathAnnot=${pathStringTie}
    export prefixAnnot=FilteredTranscripts_StringTie_Ensembl${release}
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

perl ${pathScripts}/extract.nonoverlapping.regions.pl --pathAnnotGTF=${pathAnnot}/${prefixAnnot}.gtf --pathRepeats=NA --pathRetrogenes=NA --pathOutput=${pathResults}/ExonCoordinates_ExcludingOverlapOtherGenes_${prefixAnnot}.txt

##############################################################

perl ${pathScripts}/extract.nonoverlapping.regions.pl --pathAnnotGTF=${pathAnnot}/${prefixAnnot}.gtf --pathRepeats=${pathRepeatMasker}/RepeatMasker_UCSC.txt.gz --pathRetrogenes=${pathRetrogenes} --pathOutput=${pathResults}/ExonCoordinates_ExcludingOverlapOtherGenesRepeatsRetrogenes_${prefixAnnot}.txt

##############################################################
