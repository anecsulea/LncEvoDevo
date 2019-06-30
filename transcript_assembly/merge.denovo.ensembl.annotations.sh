#!/bin/bash

export species=$1

#############################################################################

export path=LncEvoDevo

export pathDocs=${path}/docs
export pathAnnot=${path}/data/ensembl_annotations/${species}
export pathResults=${path}/results/stringtie_assembly/${species}/combined
export pathScripts=${path}/scripts/transcript_assembly

export release=94

#############################################################################

perl ${pathScripts}/merge.denovo.ensembl.annotations.pl --pathAssembledGTF=${pathResults}/assembled_transcripts.gtf --pathEnsemblGTF=${pathAnnot}/FilteredTranscripts_Ensembl${release}.gtf --pathSelectedTranscripts=${pathResults}/selected_transcripts_denovo_vs_Ensembl${release}.txt --pathOutputGTF=${pathResults}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathOutputLog=${pathResults}/log_merge_denovo_Ensembl${release}.txt

#############################################################################
