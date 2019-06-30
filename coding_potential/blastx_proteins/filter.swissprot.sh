#!/bin/bash

##############################################################

export path=LncEvoDevo
export pathSwissProt=${path}/data/SwissProt
export pathScripts=${path}/scripts/coding_potential/blastx_proteins

##############################################################

perl ${pathScripts}/filter.swissprot.pl --pathSwissProt=${pathSwissProt}/uniprot_sprot.fasta.gz --acceptedScores=1,2,3 --pathOutput=${pathSwissProt}/uniprot_sprot_filtered.fasta 

makeblastdb -in ${pathSwissProt}/uniprot_sprot_filtered.fasta -out ${pathSwissProt}/uniprot_sprot_filtered -dbtype prot -max_file_sz 3000000000B

##############################################################
