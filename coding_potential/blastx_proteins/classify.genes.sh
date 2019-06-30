#!/bin/bash

export sp=$1
export annot=$2
export db=$3

#####################################################################

export path=LncEvoDevo

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathRepeats=${path}/results/overlap_repeats/${sp}
export pathResults=${path}/results/coding_potential/${sp}
export pathScripts=${path}/scripts/coding_potential/blastx_proteins

#####################################################################

export release=94

#####################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathAnnot=${pathEnsembl}
    export suffix=FilteredTranscripts_Ensembl${release}
fi

if [ ${annot} = "StringTie" ]; then
    export pathAnnot=${pathStringTie}
    export suffix=FilteredTranscripts_StringTie_Ensembl${release}
fi


##############################################################

## protein overlap: 50 aa
## no minimum length fraction: some proteins can be short compared to UTRs

perl ${pathScripts}/classify.genes.pl --pathGTF=${pathAnnot}/${suffix}.gtf --pathFastacDNAs=${pathAnnot}/${suffix}_rm.fa --pathBlastX=${pathResults}/${suffix}_rm_vs_${db}.blastx.out --maxEValue=0.001 --minPCIdentity=40 --minFractionOverlap=0 --minLengthOverlap=150  --pathOutputGenes=${pathResults}/GeneClassification_${suffix}_rm_vs_${db}.txt  --pathOutputTranscripts=${pathResults}/TranscriptClassification_BlastX_${suffix}_rm_vs_${db}.txt 

#####################################################################
