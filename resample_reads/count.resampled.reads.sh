#!/bin/bash

export species=$1
export annot=$2

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

export release=94

###############################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathGTF=${pathAnnot}/FilteredTranscripts_Ensembl${release}.gtf
    export suffixOutput=FilteredTranscripts_Ensembl${release}.txt
fi

if [ ${annot} = "StringTie" ]; then
    export pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf
    export suffixOutput=FilteredTranscripts_StringTie_Ensembl${release}.txt
fi


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
        
    if [ -e ${pathResults}/${sample} ]; then
	echo "path output already there"
    else
	mkdir  ${pathResults}/${sample}
    fi
        
    if [ -e ${pathHisat}/${sample}/resampled_reads.sam ]; then
	## check if already done
	
	if [ -e ${pathResults}/${sample}/ReadCounts_${annot}.txt ]; then
	    echo "already done"
	else	
	    
	    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_count
	    
	  
	    
	    echo "Rscript --vanilla ${pathScripts}/count.resampled.reads.subread.R in2p3 ${annot} ${release} ${species} ${sample} 1"  >>  ${pathScripts}/bsub_script_count
	    
	    # echo "rm -r ${pathHisat}/${sample}/resampled_reads.sam" >> ${pathScripts}/bsub_script_count
	    
	  
	    qsub -q long -l s_rss=4G,sps=1  -o ${pathScripts}/std_output_count_resampled_${species}.txt -e ${pathScripts}/std_error_count_resampled_${species}.txt ${pathScripts}/bsub_script_count
	    
	fi
done 

###############################################################################
