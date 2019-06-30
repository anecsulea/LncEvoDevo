#!/bin/bash

export ref=$1
export tg=$2
export annot=$3

####################################################################################

export path=LncEvoDevo

####################################################################################

export pathLiftOver=${path}/data/genome_alignments
export pathProjections=${path}/results/exon_projections
export pathUCSC=${path}/data/UCSC_sequences
export pathEnsembl=${path}/data/ensembl_annotations
export pathStringTie=${path}/results/stringtie_assembly
export pathScripts=${path}/scripts/ortho_genes/exon_projections

export release=94

####################################################################################

if [ ${annot} = "Ensembl" ]; then
    export prefix=ExonBlocks_FilteredTranscripts_Ensembl${release}
    export pathExons=${pathEnsembl}/${ref}
fi

if [ ${annot} = "StringTie" ]; then
    export prefix=ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}
    export pathExons=${pathStringTie}/${ref}/combined
fi

####################################################################################

if [ ${ref} = "Chicken" ]||[ ${tg} = "Chicken" ]; then
    export minsizeratio=0.33
    export maxsizeratio=3
else
    export minsizeratio=0.5
    export maxsizeratio=2
fi

####################################################################################

if [ ${ref} = ${tg} ]; then
    echo "cannot project from "${ref}" to "${tg}
    exit
fi

####################################################################################

perl ${pathScripts}/filter.projected.exon.blocks.pl --pathExonBlocks=${pathExons}/${prefix}.txt --pathChromosomeCorrespondence=${pathUCSC}/${tg}/chromosomes_Ensembl_UCSC.txt --pathProjectedExons=${pathProjections}/From${ref}_To${tg}_${prefix}.bed --minSizeRatio=${minsizeratio} --maxSizeRatio=${maxsizeratio} --pathOutputFilteredExons=${pathProjections}/From${ref}_To${tg}_${prefix}_FilteredProjectedExons_Step1.txt  --pathOutputRejectedExons=${pathProjections}/From${ref}_To${tg}_${prefix}_RejectedProjectedExons_Step1.txt --pathOutputLog=${pathScripts}/logs/log_filter_exon_blocks_from${ref}_to${tg}.txt

####################################################################################
