#!/bin/bash

export sp=$1
export annot=$2

##############################################################

export path=LncEvoDevo

###########################################################################

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathBlastDB=${path}/data/blastn_databases/${sp}
export pathGenome=${path}/data/genome_indexes/${sp}
export pathScripts=${path}/scripts/ortho_genes/whole_genome_alignments

export release=94
export genomerelease=87

##############################################################

if [ ${annot} = "Ensembl" ]; then
    export pathAnnot=${pathEnsembl}
    export suffix=FilteredTranscripts_Ensembl${release}
fi

if [ ${annot} = "StringTie" ]; then
    export pathAnnot=${pathStringTie}
    export suffix=FilteredTranscripts_StringTie_Ensembl${release}
fi

##############################################################

perl ${pathScripts}/extract.cDNA.sequences.exon.blocks.pl --pathExonBlocks=${pathAnnot}/ExonBlocks_${suffix}.txt  --pathGenomeSequence=${pathGenome}/genome_ensembl${genomerelease}.fa --pathOutput=${pathAnnot}/ExonBlocks_${suffix}_cDNASequences.fa

###############################################################
