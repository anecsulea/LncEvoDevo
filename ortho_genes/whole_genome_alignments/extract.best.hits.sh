#!/bin/bash

export sp1=$1
export sp2=$2
export annot=$3
export aln=$4

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
    export prefix=FilteredTranscripts_Ensembl${release}
    export pathExons=${pathEnsembl}
    export suffix=""
fi

if [ ${annot} = "StringTie" ]; then
    export prefix=FilteredTranscripts_StringTie_Ensembl${release}
    export pathExons=${pathStringTie}
    export suffix=/combined
fi

####################################################################################

if [ ${sp1} = ${sp2} ]; then
    echo "cannot project from "${sp1}" to "${sp2}
    exit
fi

#################################################################################

perl ${pathScripts}/extract.best.hits.pl --species1=${sp1} --species2=${sp2} --pathExonBlocks1=${pathExons}/${sp1}${suffix}/ExonBlocks_${prefix}.txt --pathExonBlocks2=${pathExons}/${sp2}${suffix}/ExonBlocks_${prefix}.txt --pathUTR1=${pathExons}/${sp1}${suffix}/Putative3UTR_MaxDist1kb_${prefix}.txt --pathUTR2=${pathExons}/${sp2}${suffix}/Putative3UTR_MaxDist1kb_${prefix}.txt --pathAlignmentStats=${pathResults}/AlignmentStatistics_ExonBlocks_${aln}_${sp1}_${sp2}_ExonBlocks_${prefix}.txt  --minAlignedFraction=0 --minRatioSecondBest=1.1 --pathBestHits12=${pathResults}/BestHits_${aln}_${sp1}_${sp2}_${annot}.txt --pathBestHits21=${pathResults}/BestHits_${aln}_${sp2}_${sp1}_${annot}.txt --pathReciprocalBestHits=${pathResults}/ReciprocalBestHits_${aln}_${sp1}_${sp2}_${annot}.txt

#################################################################################
