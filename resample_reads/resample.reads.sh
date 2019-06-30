#!/bin/bash

export sp=$1

############################################################################

export path=LncEvoDevo

export pathHisat=${path}/results/hisat/${sp}
export pathScripts=${path}/scripts/resample_reads

############################################################################

for sample in `cut -f 1 ${pathHisat}/nb_resampled_unique_reads_noMT.txt | sed '1d' `
do
    if [ -e ${pathHisat}/${sample}/accepted_hits.bam ]; then

	if [ -e ${pathHisat}/${sample}/resampled_read_ids_noMT ]; then
	    echo "dir output already there"
	else
	    mkdir ${pathHisat}/${sample}/resampled_read_ids_noMT
	fi
	
	echo "#!/bin/bash" >  ${pathScripts}/bsub_script_resample
	
	
	echo "R CMD BATCH '--args sp=\"${sp}\" sample=\"${sample}\"' ${pathScripts}/resample.reads.R ${pathScripts}/resample.reads.${sp}.${sample}.Rout">>  ${pathScripts}/bsub_script_resample
	
    fi
done

############################################################################



