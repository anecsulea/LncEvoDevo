#!/bin/bash

export sp=$1

############################################################################

export path=LncEvoDevo

export pathHisat=${path}/results/hisat/${sp}
export pathScripts=${path}/scripts/resample_reads

############################################################################

for sample in `ls ${pathHisat} | grep -v txt `
do
    # if [ -e ${pathHisat}/${sample}/unique_read_ids_noMT.txt.gz ]; then
    if [ -e ${pathHisat}/${sample}/unique_read_ids.txt.gz ]||[ -e ${pathHisat}/${sample}/unique_read_ids.txt ]; then
	echo ${sample} "already done"
    else
	if [ -e ${pathHisat}/${sample}/accepted_hits.bam ]; then
	    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_uniqueid
	   
	    
	    echo "perl ${pathScripts}/extract.unique.read.ids.pl --pathAlignments=${pathHisat}/${sample}/accepted_hits.bam --forbiddenChromo=MT --pathOutput=${pathHisat}/${sample}/unique_read_ids_noMT.txt --pathOutputCounts=${pathHisat}/${sample}/nb_unique_reads_noMT.txt">>  ${pathScripts}/bsub_script_uniqueid
	    
	    echo "gzip ${pathHisat}/${sample}/unique_read_ids_noMT.txt">>  ${pathScripts}/bsub_script_uniqueid
	    
	 
	    qsub -q long -l s_rss=4G,sps=1 -o ${pathScripts}/std_output_uniqueid_${sp}.txt -e ${pathScripts}/std_error_uniqueid_${sp}.txt ${pathScripts}/bsub_script_uniqueid
		
	fi
    fi
done

############################################################################



