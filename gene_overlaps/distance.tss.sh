#!/bin/bash

export sp=$1

#########################################################################

export path=LncEvoDevo

export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathResults=${path}/results/gene_overlaps/${sp}
export pathScripts=${path}/scripts/gene_overlaps

export release=94

#########################################################################

perl ${pathScripts}/distance.tss.pl --pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --maxDistance=1000 --pathOutput=${pathResults}/DistanceAntisenseTSS_MaxDist1kb_FilteredTranscripts_StringTie_Ensembl${release}.txt


perl ${pathScripts}/distance.tss.pl --pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --maxDistance=5000 --pathOutput=${pathResults}/DistanceAntisenseTSS_MaxDist5kb_FilteredTranscripts_StringTie_Ensembl${release}.txt


#########################################################################
