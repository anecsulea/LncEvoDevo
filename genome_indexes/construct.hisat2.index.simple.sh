#!/bin/bash 

export sp=$1

#####################################################################

export path=LncEvoDevo
export pathHisat2=${path}/data/genome_indexes/${sp}
export pathScripts=${path}/scripts/genome_indexes

## hisat 2.0.5

export release=87

######################################################################

if [ -e ${pathHisat2} ]; then
    echo "path output already there"
else
    mkdir ${pathHisat2}
fi

######################################################################

echo "#!/bin/bash" >  ${pathScripts}/bsub_script_hisat2
echo "#PBS -o ${pathScripts}/std_output_hisat2_${sp}.txt" >>  ${pathScripts}/bsub_script_hisat2
echo "#PBS -e ${pathScripts}/std_error_hisat2_${sp}.txt" >>  ${pathScripts}/bsub_script_hisat2

echo "hisat2-build -p 4 --seed 19 ${pathHisat2}/genome_ensembl${release}.fa ${pathHisat2}/genome_ensembl${release}"  >> ${pathScripts}/bsub_script_hisat2 

qsub -q q1day -l nodes=1:ppn=4,mem=10gb ${pathScripts}/bsub_script_hisat2

######################################################################
