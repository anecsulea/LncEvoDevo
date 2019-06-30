#!/bin/bash

export sp=$1

#########################################################################

export path=LncEvoDevo


export pathCGI=${path}/data/CpGIslands/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathResults=${path}/results/gene_overlaps/${sp}
export pathScripts=${path}/scripts/gene_overlaps

export release=94

#########################################################################

perl ${pathScripts}/overlap.tss.cgi.pl --pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathCpGIslands=${pathCGI}/CpGIslands_UCSC.txt --maxDistance=5000 --pathOutput=${pathResults}/OverlapCGI_MaxDist5kb_FilteredTranscripts_StringTie_Ensembl${release}.txt


perl ${pathScripts}/overlap.tss.cgi.pl --pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathCpGIslands=${pathCGI}/CpGIslands_UCSC.txt --maxDistance=1000 --pathOutput=${pathResults}/OverlapCGI_MaxDist1kb_FilteredTranscripts_StringTie_Ensembl${release}.txt



perl ${pathScripts}/overlap.tss.cgi.pl --pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathCpGIslands=${pathCGI}/CpGIslands_UCSC.txt --maxDistance=500 --pathOutput=${pathResults}/OverlapCGI_MaxDist500bp_FilteredTranscripts_StringTie_Ensembl${release}.txt



perl ${pathScripts}/overlap.tss.cgi.pl --pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathCpGIslands=${pathCGI}/CpGIslands_UCSC.txt --maxDistance=0 --pathOutput=${pathResults}/OverlapCGI_StrictOverlap_FilteredTranscripts_StringTie_Ensembl${release}.txt

#########################################################################

