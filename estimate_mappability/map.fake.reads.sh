#/bin/bash

export species=$1

##############################################################

export path=LncEvoDevo
export pathResults=${path}/results/mappability/${species}
export pathGenomeIndexes=${path}/data/genome_indexes/${species}
export pathScripts=${path}/scripts/estimate_mappability

export release=87
export pathIndex=${pathGenomeIndexes}/genome_ensembl${release}

##############################################################

for chr in {1..32} X Y Z W
do
    
    if [ -e ${pathResults}/fake_reads_chr${chr}.sam.gz ]; then
	echo "already done"
    else
	if [ -e ${pathResults}/fake_reads_posinfo_chr${chr}.fa.gz ]; then
	    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_map
	    echo "#PBS -o ${pathScripts}/std_output_map_${species}_${chr}.txt" >>  ${pathScripts}/bsub_script_map
	    echo "#PBS -e ${pathScripts}/std_error_map_${species}_${chr}.txt" >>  ${pathScripts}/bsub_script_map
	    echo "source /panhome/necsulea/.bashrc" >>  ${pathScripts}/bsub_script_map
	    
	    export pathLocal="./map_fakereads_${chr}_${species}"
	    
	    echo "if [ -e ${pathLocal} ]; then" >> ${pathScripts}/bsub_script_map
	    echo "echo 'path exists'">> ${pathScripts}/bsub_script_map
	    echo "else">> ${pathScripts}/bsub_script_map
	    echo "mkdir ${pathLocal}">> ${pathScripts}/bsub_script_map
	    echo "fi" >> ${pathScripts}/bsub_script_map
	    
	    echo "hisat2 -f --seed 19 -p 12 -x ${pathIndex} -U ${pathResults}/fake_reads_posinfo_chr${chr}.fa.gz -S ${pathLocal}/fake_reads_chr${chr}.sam  --max-intronlen 1000000 --dta --no-unal " >> ${pathScripts}/bsub_script_map
	
	    echo "gzip ${pathLocal}/fake_reads_chr${chr}.sam" >> ${pathScripts}/bsub_script_map
	    echo "mv ${pathLocal}/fake_reads_chr${chr}.sam.gz ${pathResults}/" >> ${pathScripts}/bsub_script_map
	    
	    echo "rm -r ${pathLocal}" >> ${pathScripts}/bsub_script_map
	    
	    qsub -q q1day -l nodes=1:ppn=12,mem=5gb ${pathScripts}/bsub_script_map
	fi
    fi
done

##############################################################


