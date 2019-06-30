#!/bin/bash

export ref=$1
export tg=$2
export annot=$3

######################################################################

export path=LncEvoDevo

######################################################################

export pathEnsembl=${path}/data/ensembl_annotations
export pathStringTie=${path}/results/stringtie_assembly
export pathGenomeAlignments=${path}/data/genome_alignments
export pathResults=${path}/results/exon_projections
export pathScripts=${path}/scripts/exon_projections

export release=94

#####################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathExons=${pathEnsembl}/${ref}/bed_parts
    export prefixExons=ExonBlocks_FilteredTranscripts_Ensembl${release}
fi

if [ ${annot} = "StringTie" ]; then
    export pathExons=${pathStringTie}/${ref}/combined/bed_parts
    export prefixExons=ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}
fi

######################################################################

if [ -e ${pathResults}/From${ref}_To${tg}_${prefixExons}.bed ]; then
    echo "already done"
else
    for file in `ls ${pathExons} | grep ${prefixExons} | grep part`
    do
	export prefix=`basename ${file} .bed`
	echo ${file} ${prefix}
	
	if [ -e ${pathResults}/From${ref}_To${tg}_${prefix}.bed ]; then
	    cat ${pathResults}/From${ref}_To${tg}_${prefix}.bed >> ${pathResults}/From${ref}_To${tg}_${prefixExons}.bed
	    cat ${pathResults}/From${ref}_To${tg}_${prefix}.unmapped >> ${pathResults}/From${ref}_To${tg}_${prefixExons}.unmapped
	else
	    echo "cannot find liftOver results"
	fi
    done
fi

#####################################################################
