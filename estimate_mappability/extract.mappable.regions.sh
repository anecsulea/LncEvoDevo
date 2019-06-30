#/bin/bash

export species=$1

##############################################################

export path=LncEvoDevo
export pathResults=${path}/results/mappability/${species}
export pathScripts=${path}/scripts/estimate_mappability

##############################################################

for chr in {1..32} X Y Z W 
do
    
    if [ -e ${pathResults}/fake_reads_chr${chr}_mappedregions.txt ]; then
	echo "already done"
    else
	if [ -e ${pathResults}/fake_reads_posinfo_chr${chr}.fa.gz ]&&[ -e ${pathResults}/fake_reads_chr${chr}.sam.gz ]; then
    	   
	    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_estmap
	    echo "#PBS -o ${pathScripts}/std_output_estmap_${species}_${chr}.txt" >>  ${pathScripts}/bsub_script_estmap
	    echo "#PBS -e ${pathScripts}/std_error_estmap_${species}_${chr}.txt" >>  ${pathScripts}/bsub_script_estmap
	    echo "source /panhome/necsulea/.bashrc" >>  ${pathScripts}/bsub_script_estmap
	    
	    echo "perl ${pathScripts}/extract.mappable.regions.pl --pathReads=${pathResults}/fake_reads_posinfo_chr${chr}.fa.gz --pathAlignment=${pathResults}/fake_reads_chr${chr}.sam.gz --readLength=101 --pathOutput=${pathResults}/fake_reads_chr${chr}_mappedregions.txt" >> bsub_script_estmap
	    
	    qsub -q q1day -l nodes=1:ppn=1,mem=25gb ${pathScripts}/bsub_script_estmap
	    
	fi
    fi

done


##############################################################
