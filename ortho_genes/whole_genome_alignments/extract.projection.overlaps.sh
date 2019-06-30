#!/bin/bash

export sp1=$1
export sp2=$2
export annot=$3

####################################################################################

export path=LncEvoDevo

####################################################################################

export pathLiftOver=${path}/data/genome_alignments
export pathProjections=${path}/results/exon_projections
export pathStringTie=${path}/results/stringtie_assembly
export pathEnsembl=${path}/data/ensembl_annotations
export pathUCSC=${path}/data/UCSC_sequences
export pathResults=${path}/results/ortho_genes/whole_genome_alignments
export pathScripts=${path}/scripts/ortho_genes/whole_genome_alignments

export release=94

####################################################################################

if [ ${annot} = "Ensembl" ]; then
    export prefix=ExonBlocks_FilteredTranscripts_Ensembl${release}
    export pathExons=${pathEnsembl}
    export suffix=""
fi

if [ ${annot} = "StringTie" ]; then
    export prefix=ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}
    export pathExons=${pathStringTie}
    export suffix="/combined"
fi

####################################################################################

if [ ${sp1} = ${sp2} ]; then
    echo "cannot project from "${sp1}" to "${sp2}
    exit
fi

####################################################################################

perl ${pathScripts}/extract.projection.overlaps.pl --species1=${sp1} --species2=${sp2} --pathExonBlocks1=${pathExons}/${sp1}${suffix}/${prefix}.txt --pathExonBlocks2=${pathExons}/${sp2}${suffix}/${prefix}.txt --pathProjectedExons12=${pathProjections}/From${sp1}_To${sp2}_${prefix}_FilteredProjectedExons_Step2.txt --pathProjectedExons21=${pathProjections}/From${sp2}_To${sp1}_${prefix}_FilteredProjectedExons_Step2.txt --pathProjectionMap12=${pathResults}/ProjectionMap_From${sp1}_To${sp2}_${prefix}.txt --pathProjectionMap21=${pathResults}/ProjectionMap_From${sp2}_To${sp1}_${prefix}.txt --pathGeneClusters=${pathResults}/ProjectionClusters_${sp1}_${sp2}_${prefix}.txt 

####################################################################################
