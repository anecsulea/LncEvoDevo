#/bin/bash

export species=$1

##############################################################

export path=LncEvoDevo
export pathSequence=${path}/data/genome_indexes/${species}
export pathResults=${path}/results/mappability/${species}
export pathScripts=${path}/scripts/estimate_mappability

export release=87

##############################################################

if [ -e ${pathSequence}/genome_ensembl${release}.fa ]; then
    export suffix=fa
else
    if [ -e ${pathSequence}/genome_ensembl${release}.fa.gz ]; then
	export suffix=fa.gz
    else
	echo "cannot find fasta file"
	exit
    fi
fi

##############################################################

for chr in {1..32} X Y Z W
do
    
    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_estmap
    echo "#PBS -o ${pathScripts}/std_output_estmap.txt" >>  ${pathScripts}/bsub_script_estmap
    echo "#PBS -e ${pathScripts}/std_error_estmap.txt" >>  ${pathScripts}/bsub_script_estmap
    echo "source /panhome/necsulea/.bashrc" >>  ${pathScripts}/bsub_script_estmap

    export pathLocal="./fakereads_${chr}_${species}"

    echo "if [ -e ${pathLocal} ]; then" >> ${pathScripts}/bsub_script_estmap
    echo "echo 'path exists'">> ${pathScripts}/bsub_script_estmap
    echo "else">> ${pathScripts}/bsub_script_estmap
    echo "mkdir ${pathLocal}">> ${pathScripts}/bsub_script_estmap
    echo "fi" >> ${pathScripts}/bsub_script_estmap
    
    echo "perl ${pathScripts}/generate.fake.reads.pl --pathGenome=${pathSequence}/genome_ensembl${release}.${suffix} --readLength=101 --step=5 --chr=${chr} --pathOutput=${pathLocal}/fake_reads_posinfo_chr${chr}.fa" >> ${pathScripts}/bsub_script_estmap

    echo "gzip ${pathLocal}/fake_reads_posinfo_chr${chr}.fa" >> ${pathScripts}/bsub_script_estmap
    
    echo "mv ${pathLocal}/fake_reads_posinfo_chr${chr}.fa.gz ${pathResults}/">> ${pathScripts}/bsub_script_estmap

    echo "rm -r ${pathLocal}" >> ${pathScripts}/bsub_script_estmap

    qsub -q q1hour -l nodes=1:ppn=1,mem=5gb ${pathScripts}/bsub_script_estmap
    
done

##############################################################



