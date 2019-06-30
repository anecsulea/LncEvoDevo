#!/bin/bash

export sp=$1

############################################################################

export path=LncEvoDevo

export pathHisat=${path}/results/hisat/${sp}
export pathScripts=${path}/scripts/resample_reads

##############################################################################

if [ -e ${pathHisat}/nb_unique_reads_noMT.txt ]; then
    echo "cleaning up previous file"
    rm ${pathHisat}/nb_unique_reads_noMT.txt
fi

##############################################################################

for sample in `ls ${pathHisat} | grep -v txt`
do
    if [ -e ${pathHisat}/${sample}/accepted_hits.bam ]; then
	echo -n ${sample}$'\t' >> ${pathHisat}/nb_unique_reads_noMT.txt
	cat  ${pathHisat}/${sample}/nb_unique_reads_noMT.txt  >> ${pathHisat}/nb_unique_reads_noMT.txt
    fi
done

##############################################################################
