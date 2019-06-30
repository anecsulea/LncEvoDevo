#!/bin/bash

export sp=$1
export alntype=$2
export windowSize=$3

##################################################################

export path=LncEvoDevo
export pathMatrices=${path}/results/CSF/matrices/${alntype}/${sp}
export pathResults=${path}/results/CSF/genome_scores/${alntype}/${sp}
export pathScripts=${path}/scripts/coding_potential/CSF
export pathGenomeAlignments=${path}/data/genome_alignments/${alntype}
export pathGeneticCode=${pathScripts}/standard_genetic_code.txt 

##################################################################

if [ ${sp} = "Mouse" ]; then
   export refsp="mm10"
fi

if [ ${sp} = "Chicken" ]; then
   export refsp="galGal4"
fi

if [ ${sp} = "Rat" ]; then
   export refsp="rn6"
fi

if [ ${sp} = "Human" ]; then
   export refsp="hg38"
fi

##################################################################

export penalty=0.001
export pseudofreq=0.00001

export maxsize=300000000

if [ ${sp} = "Opossum" ]; then
    export maxsize=800000000
fi 

echo ${refsp}" max size "${maxsize}

##################################################################

for chr in {1..32} X Y Z W
do
    if [ -e ${pathGenomeAlignments}/chr${chr}.maf.gz ]; then
	
	export start="0"
	
	while [ ${start} -lt ${maxsize} ]
	do
	    export end=$[${start}+20000000]

	    echo ${chr} ${start} ${end}
	    
	    if [ -e ${pathResults}/CSF_PositiveScores_windowSize${windowSize}bp_penalty${penalty}_pseudofreq${pseudofreq}_chr${chr}_start${start}_end${end}.txt ]&&[ -e ${pathResults}/CSF_CoveredRegions_windowSize${windowSize}bp_penalty${penalty}_pseudofreq${pseudofreq}_chr${chr}_start${start}_end${end}.txt ]; then
		echo "already done"
	    else
		
		echo "#!/bin/bash" >   ${pathScripts}/bsub_script_scores
		echo "#PBS -o ${pathScripts}/std_output_scores_${refsp}_${chr}.txt" >>  ${pathScripts}/bsub_script_scores
		echo "#PBS -e ${pathScripts}/std_error_scores_${refsp}_${chr}.txt" >>  ${pathScripts}/bsub_script_scores
		echo "source /panhome/necsulea/.bashrc" >>   ${pathScripts}/bsub_script_scores
		
		echo "perl ${pathScripts}/csf.score.genome.pl --pathAlignment=${pathGenomeAlignments}/chr${chr}.maf.gz --refsp=${refsp} --pathSelectedInformants=${pathMatrices}/selected_informant_species.txt  --pathCodingMatrices=${pathMatrices}/matrices_CDS.txt --pathNoncodingMatrices=${pathMatrices}/matrices_introns.txt --pathGeneticCode=${pathGeneticCode} --windowSize=${windowSize}  --minInformantSpecies=3 --minInformativeCodons=5 --penalty=${penalty} --pseudofreq=${pseudofreq} --includeSymmetric=yes --start=${start} --end=${end} --pathOutputPositiveScores=${pathResults}/CSF_PositiveScores_windowSize${windowSize}bp_penalty${penalty}_pseudofreq${pseudofreq}_chr${chr}_start${start}_end${end}.txt --pathOutputCoveredRegions=${pathResults}/CSF_CoveredRegions_windowSize${windowSize}bp_penalty${penalty}_pseudofreq${pseudofreq}_chr${chr}_start${start}_end${end}.txt">> ${pathScripts}/bsub_script_scores
		
		qsub -q q1hour -l nodes=1,mem=10gb ${pathScripts}/bsub_script_scores
	    fi
	    
	    export start=$[${start}+20000000]
	done
    fi
done

##################################################################
