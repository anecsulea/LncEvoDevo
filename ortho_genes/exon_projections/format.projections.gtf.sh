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
export pathScripts=${path}/scripts/exon_projections

export release=89

####################################################################################

if [ ${annot} = "Ensembl" ]; then
    export prefix=FilteredTranscripts_Ensembl${release}
    export pathExons=${pathEnsembl}/${ref}
fi

if [ ${annot} = "StringTie" ]; then
    export prefix=FilteredTranscripts_StringTie_Ensembl${release}
    export pathExons=${pathStringTie}/${ref}/combined
fi

####################################################################################

if [ -e ${pathProjections}/From${ref}_To${tg}_${prefix}_FilteredProjectedExons_Step2.txt ]; then
    perl ${pathScripts}/format.projections.gtf.pl --pathProjectedExons=${pathProjections}/From${ref}_To${tg}_${prefix}_FilteredProjectedExons_Step2.txt --pathOutput=${pathProjections}/From${ref}_To${tg}_${prefix}_FilteredProjectedExons_Step2.gtf
fi

####################################################################################
