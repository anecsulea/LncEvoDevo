#!/bin/bash

##############################################################

export path=LncEvoDevo
export pathPfam=${path}/data/Pfam
export pathSwissProt=${path}/data/SwissProt
export pathScripts=${path}/scripts/coding_potential/blastx_proteins

##############################################################

makeblastdb -in ${pathPfam}/Pfam-A.fasta -out ${pathPfam}/Pfam-A -dbtype prot -max_file_sz 10GB


#makeblastdb -in ${pathSwissProt}/uniprot_sprot_filtered.fasta -out ${pathSwissProt}/uniprot_sprot_filtered -dbtype prot -max_file_sz 10GB


##############################################################
