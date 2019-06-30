#!/bin/bash

export type=$1

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

####################################################################################

export pathRBH=""

for file in `ls ${pathResults} | grep ReciprocalBestHits | grep ${type}`
do
    export pathRBH=${pathResults}/${file},${pathRBH}
done

echo ${pathRBH}

####################################################################################

perl ${pathScripts}/extract.gene.families.pl --pathsReciprocalBestHits=${pathRBH} --pathOutput1to1Families=${pathResults}/GeneFamilies_1to1_${type}.txt --pathOutputMultipleFamilies=${pathResults}/GeneFamilies_Multiple_${type}.txt 

####################################################################################
