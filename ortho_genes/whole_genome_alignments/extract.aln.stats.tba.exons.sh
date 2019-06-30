#!/bin/bash

export sp1=$1
export sp2=$2
export annot=$3

####################################################################

export path=LncEvoDevo

export pathStringTie=${path}/results/stringtie_assembly
export pathGenomeAlignments=${path}/data/genome_alignments
export pathGenomeSequences=${path}/data/genome_indexes
export pathResults=${path}/results/ortho_genes/whole_genome_alignments
export pathScripts=${path}/scripts/ortho_genes/whole_genome_alignments

export release=94

####################################################################

if [ ${annot} = "StringTie" ]; then
    export prefixExons=ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}
fi


if [ ${annot} = "Ensembl" ]; then
    export prefixExons=ExonBlocks_FilteredTranscripts_Ensembl${release}
fi

####################################################################
        
perl ${pathScripts}/extract.aln.stats.tba.exons.pl --species1=${sp1} --species2=${sp2} --pathClusters=${pathResults}/ProjectionClusters_${sp1}_${sp2}_ExonBlocks_FilteredTranscripts_StringTie_Ensembl94.txt --dirTBA=${pathResults}/tba_alignments_projection_clusters/${sp1}_${sp2}_StringTie/ --minAlignmentLength=0 --pathOutput=${pathResults}/AlignmentStatistics_ExonBlocks_TBA_${sp1}_${sp2}_${prefixExons}.txt 

####################################################################
