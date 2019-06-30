#!/bin/bash

export ref=$1
export tg=$2
export annot=$3

####################################################################################

export path=LncEvoDevo

####################################################################################

export pathLiftOver=${path}/data/genome_alignments
export pathProjections=${path}/results/exon_projections
export pathStringTie=${path}/results/stringtie_assembly
export pathUCSC=${path}/data/UCSC_sequences
export pathEnsembl=${path}/data/ensembl_annotations
export pathSynteny=${path}/results/ortho_genes/synteny
export pathScripts=${path}/scripts/ortho_genes/exon_projections

export release=94

####################################################################################

if [ ${annot} = "Ensembl" ]; then
    export prefix=ExonBlocks_FilteredTranscripts_Ensembl${release}
    export suffixAnnot=FilteredTranscripts_Ensembl${release}
    export pathExons=${pathEnsembl}/${ref}
fi

if [ ${annot} = "StringTie" ]; then
    export prefix=ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}
    export suffixAnnot=FilteredTranscripts_StringTie_Ensembl${release}
    export pathExons=${pathStringTie}/${ref}/combined
fi

####################################################################################

if [ ${ref} = ${tg} ]; then
    echo "cannot project from "${ref}" to "${tg}
    exit
fi

####################################################################################

perl ${pathScripts}/filter.projected.genes.exon.blocks.pl --pathExonBlocks=${pathExons}/${prefix}.txt --pathProjectedExons=${pathProjections}/From${ref}_To${tg}_${prefix}_FilteredProjectedExons_Step1.txt --maxIntronSizeRatio=100 --maxAddedIntronSize=1000000 --pathSyntenyPredictions=NA --syntenyRange=1000000 --pathOutputLog=${pathScripts}/logs/log_filter_genes_exon_blocks_from${ref}_to${tg}.txt --pathOutputFilteredExons=${pathProjections}/From${ref}_To${tg}_${prefix}_FilteredProjectedExons_Step2.txt --pathOutputRejectedExons=${pathProjections}/From${ref}_To${tg}_${prefix}_RejectedProjectedExons_Step2.txt 

####################################################################################
