#!/bin/bash

export sp=$1

#####################################################################

export path=LncEvoDevo

export pathAnnot=${path}/results/stringtie_assembly/${sp}/combined
export pathScripts=${path}/scripts/transcript_assembly

export release=94

#####################################################################

for suffix in assembled_transcripts FilteredTranscripts_StringTie_Ensembl${release}
do
    perl ${pathScripts}/make.exon.blocks.gtf.pl --pathGTF=${pathAnnot}/${suffix}.gtf --collapseDistance=0  --pathOutputExonBlocks=${pathAnnot}/ExonBlocks_${suffix}.txt
done

#################################################################################

    
