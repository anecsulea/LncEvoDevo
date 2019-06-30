#!/bin/bash

export species=$1
export type=$2

####################################################################################

export path=LncEvoDevo
    

####################################################################################

export pathDocs=${path}/docs
export pathRNASeq=${path}/data/RNASeq/${species}
export pathResults=${path}/results/expression_estimation/${species}
export pathIndexes=${path}/results/kallisto_indexes/${species}
export pathScripts=${path}/scripts/expression_estimation

export release=94

####################################################################################

for file in `ls ${pathRNASeq} | grep _cutadapt.fastq.gz`
do
    export sample=`basename ${file} _trimmed_cutadapt.fastq.gz`

    if [ -e ${pathResults}/${sample} ]; then 
	echo "path exists"
    else
	mkdir ${pathResults}/${sample}
    fi
    
    if [ ${species} = "Chicken" ]; then
	export library=`grep ^${sample}$'\t' ${pathDocs}/LibraryType_Chicken.txt | cut -f 2`
    else
	if [ ${sample} = "Sertoli" ]||[ ${sample} = "Spermatogonia" ]||[ ${sample} = "Spermatocytes" ]||[ ${sample} = "Spermatids" ]||[ ${sample} = "Spermatozoa" ]; then
	    export library=fr-secondstrand
	else
	    export library=fr-firststrand
	fi
    fi
    
    export ss=""
    
    if [ ${library} = "fr-secondstrand" ]; then
	export ss="--fr-stranded"
    fi
    
    if [ ${library} = "fr-firststrand" ]; then
	export ss="--rf-stranded"
    fi
    
	
    echo ${sample} ${library} ${ss}
    
    if [ -e ${pathResults}/${sample}/kallisto_${type}/abundance.tsv ]; then
	echo "already done"
    else
	if [ -e ${pathResults}/${sample}/kallisto_${type} ]; then
	    echo "output dir already there"
	else
	    mkdir ${pathResults}/${sample}/kallisto_${type}
	fi
	
	echo "#!/bin/bash" >  ${pathScripts}/bsub_script

	echo "kallisto quant --single -l 200.0 -s 20 --bias -t 4 ${ss} -o ${pathResults}/${sample}/kallisto_${type} --index ${pathIndexes}/${type} ${pathRNASeq}/${sample}_trimmed_cutadapt.fastq.gz " >> ${pathScripts}/bsub_script
	    	

	qsub -q mc_huge -l s_rss=8G,sps=1 -pe multicores 4 -o ${pathScripts}/std_output_kallisto_${sp}.txt -e ${pathScripts}/std_error_kallisto_${sp}.txt ${pathScripts}/bsub_script

    fi
done

####################################################################################
