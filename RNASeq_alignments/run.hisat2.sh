#!/bin/bash

export species=$1
export tissue=$2

#############################################################################

export path=/pandata/necsulea/LncEvoDevo
export pathRNASeq=${path}/data/RNASeq/${species}
export pathEnsembl=${path}/data/ensembl_annotations/${species}
export pathDocs=${path}/docs
export pathResults=${path}/results/hisat/${species}
export pathGenomeIndexes=${path}/data/genome_indexes/${species}

export release=87

export pathIndex=${pathGenomeIndexes}/genome_ensembl${release}

#############################################################################

for file in `ls ${pathRNASeq} | grep fastq.gz | grep ${tissue} `
do
    export sample=`basename ${file} _trimmed_cutadapt.fastq.gz`
   
    if [ -e ${pathRNASeq}/${sample}_trimmed_cutadapt.fastq.gz ]; then
	export suffix=_trimmed_cutadapt.fastq.gz
    else
	if [ -e ${pathRNASeq}/${sample}.fastq.gz ]; then
	    export suffix=.fastq.gz
	else
	    echo "Weird! cannot find data for "${sample}
	fi
    fi
    
    if [ -e ${pathResults}/${sample}/accepted_hits.sam.gz ]||[ -e ${pathResults}/${sample}/accepted_hits.bam ]; then
    	echo "already done"
    else

	if [ -e ${pathResults}/${sample} ]; then
    	    echo "dir output already there"
	else
    	    mkdir ${pathResults}/${sample}
	fi
	
	export strand="--rna-strandness R"
	
	if [ ${species} = "Chicken" ]; then
	    export library=`grep ^${sample}$'\t' ${pathDocs}/LibraryType_Chicken.txt | cut -f 2`
	    
	    if [ ${library} = "fr-unstranded" ]; then
		export strand=""
	    fi

	    if [ ${library} = "fr-firststrand" ]; then
		export strand="--rna-strandness R"
	    fi

	    if [ ${library} = "fr-secondstrand" ]; then
		export strand="--rna-strandness F"
	    fi
	fi
	
	echo ${sample} ${strand}
	
	echo "#!/bin/bash" > bsub_script_hisat
	echo "#PBS -o std_output_hisat_${species}_${sample}.txt" >>  bsub_script_hisat
	echo "#PBS -e std_error_hisat_${species}_${sample}.txt" >>  bsub_script_hisat
	echo "source /panhome/necsulea/.bashrc" >>  bsub_script_hisat
    
    	export pathLocal=./hisat2_${species}_${sample}
	
    	echo "if [ -e ${pathLocal} ]; then">> bsub_script_hisat
    	echo "echo local path exists">> bsub_script_hisat
    	echo "else">> bsub_script_hisat
    	echo "mkdir ${pathLocal}">> bsub_script_hisat
    	echo "fi">> bsub_script_hisat
	
    	echo "hisat2 --seed 19 -p 6 -x ${pathIndex} -U ${pathRNASeq}/${sample}${suffix} -S ${pathLocal}/accepted_hits.sam ${strand} --known-splicesite-infile=${pathEnsembl}/SpliceSites_Ensembl87.txt --max-intronlen 1000000 --dta-cufflinks --no-unal --met-file ${pathLocal}/metrics.txt  --novel-splicesite-outfile ${pathLocal}/novel_splicesites.txt >& ${pathLocal}/align_summary.txt">> bsub_script_hisat

	echo "if [ -e ${pathLocal}/accepted_hits.sam ]; then">> bsub_script_hisat
	echo "rm ${pathRNASeq}/${sample}${suffix}" >> bsub_script_hisat
	echo "fi">> bsub_script_hisat

	echo "samtools sort -m 8G -@ 6 -o ${pathLocal}/accepted_hits.bam -O bam ${pathLocal}/accepted_hits.sam ">> bsub_script_hisat

	echo "rm ${pathLocal}/accepted_hits.sam " >> bsub_script_hisat

    	echo "mv ${pathLocal}/* ${pathResults}/${sample}/" >> bsub_script_hisat
	
	echo "rm -r ${pathLocal}" >> bsub_script_hisat 

	qsub -q q1day -l nodes=1:ppn=6,mem=10gb bsub_script_hisat
    fi
done

#############################################################################
