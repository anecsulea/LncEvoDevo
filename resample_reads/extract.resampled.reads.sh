#!/bin/bash

export species=$1

##############################################################################

export path=LncEvoDevo

#############################################################################

export pathAnnot=${path}/data/ensembl_annotations/${species}
export pathHisat=${path}/results/hisat/${species}
export pathStringTie=${path}/results/stringtie_assembly/${species}/combined
export pathResults=${path}/results/read_counts_resampled_noMT/${species}
export pathScripts=${path}/scripts/resample_reads
export pathDocs=${path}/docs

###############################################################################

## check output directory

if [ -e ${pathResults} ]; then
    echo "dir output already there"
else
    mkdir ${pathResults}
fi

###############################################################################

for sample in `cut -f 1 ${pathHisat}/nb_resampled_unique_reads_noMT.txt | sed '1d'`
do
    export nbreads=`grep ^${sample}$'\t' ${pathHisat}/nb_resampled_unique_reads_noMT.txt | cut -f 4`

    echo ${sample} ${nbreads}

    if [ -e ${pathHisat}/${sample}/resampled_reads.sam ]; then
	echo "already done"
    else
        if [ -e ${pathHisat}/${sample}/accepted_hits.bam ]; then

	    if [ -e ${pathHisat}/${sample}/resampled_read_ids_noMT/read_ids_${nbreads}.txt ]; then
		
		## check if already done
		
		if [ -e ${pathResults}/${sample}/ReadCounts_${annot}.txt ]; then
		    echo "already done"
		else	
	   	    
		    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_count
		    
		   
		    echo "perl ${pathScripts}/select.alignments.pl --pathAlignments=${pathHisat}/${sample}/accepted_hits.bam --pathAllReads=${pathHisat}/${sample}/unique_read_ids_noMT.txt.gz --pathSelectedReads=${pathHisat}/${sample}/resampled_read_ids_noMT/read_ids_${nbreads}.txt --pathOutput=${pathHisat}/${sample}/resampled_reads.sam " >>  ${pathScripts}/bsub_script_count
		    
		    qsub -q mc_highmem_long -l s_rss=20G,sps=1  -o ${pathScripts}/std_output_extract_resampled_${species}.txt -e ${pathScripts}/std_error_extract_resampled_${species}.txt ${pathScripts}/bsub_script_count
		    
	    fi
	fi
    fi
done 

###############################################################################
