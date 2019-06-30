#!/bin/bash

export sp=$1
export aln=$2

############################################################################################

export path=LncEvoDevo
export pathTrainingDataset=${path}/results/CSF/training_dataset/${sp}
export pathMatrices=${path}/results/CSF/matrices/${aln}/${sp}
export pathAln=${path}/data/genome_alignments/${aln}
export pathScripts=${path}/scripts/coding_potential/CSF

export release=84

############################################################################################

if [ ${sp} = "Mouse" ]&&[ ${aln} = "Mouse_60way" ]; then
    export refsp="mm10"
fi

if [ ${sp} = "Rat" ]&&[ ${aln} = "Human_100way" ]; then
    export refsp="rn6"
fi

if [ ${sp} = "Rat" ]&&[ ${aln} = "Rat_20way" ]; then
    export refsp="rn6"
fi

if [ ${sp} = "Chicken" ]&&[ ${aln} = "Human_100way" ]; then
    export refsp="galGal4"
fi

if [ ${sp} = "Human" ]&&[ ${aln} = "Human_100way" ]; then
    export refsp="hg38"
fi

############################################################################################
    
for chr in {1..22} X Y 
do
    if [ -e ${pathAln}/chr${chr}.maf.gz ]; then
	
	echo "#!/bin/bash" >   ${pathScripts}/bsub_script_int
	echo "#PBS -o ${pathScripts}/std_output_matrices_int_${sp}_${chr}.txt" >>  ${pathScripts}/bsub_script_int 
	echo "#PBS -e ${pathScripts}/std_error_matrices_int_${sp}_${chr}.txt" >>  ${pathScripts}/bsub_script_int 
	echo "source /panhome/necsulea/.bashrc" >>   ${pathScripts}/bsub_script_int
	
	echo "perl ${pathScripts}/compute.matrices.pl --pathAlignment=${pathAln}/chr${chr}.maf.gz --refSpecies=${refsp} --pathRegionCoordinates=${pathTrainingDataset}/SelectedIntrons.txt --minAlignedLength=30 --minAlignedSpecies=2 --pathGeneticCode=${pathScripts}/standard_genetic_code.txt --pathOutput=${pathMatrices}/matrices_introns_chr${chr}.txt" >> ${pathScripts}/bsub_script_int

	qsub -l walltime=4:00:00,nodes=1,mem=3gb ${pathScripts}/bsub_script_int

	echo "#!/bin/bash" >   ${pathScripts}/bsub_script_cds
	echo "#PBS -o ${pathScripts}/std_output_matrices_cds_${sp}_${chr}.txt" >>  ${pathScripts}/bsub_script_cds 
	echo "#PBS -e ${pathScripts}/std_error_matrices_cds_${sp}_${chr}.txt" >>  ${pathScripts}/bsub_script_cds 
	echo "source /panhome/necsulea/.bashrc" >>  ${pathScripts}/bsub_script_cds
		
	echo "perl ${pathScripts}/compute.matrices.pl --pathAlignment=${pathAln}/chr${chr}.maf.gz --refSpecies=${refsp}  --pathRegionCoordinates=${pathTrainingDataset}/SelectedCDS.txt --minAlignedLength=30 --minAlignedSpecies=2 --pathGeneticCode=${pathScripts}/standard_genetic_code.txt --pathOutput=${pathMatrices}/matrices_CDS_chr${chr}.txt" >> ${pathScripts}/bsub_script_cds

	qsub -l walltime=4:00:00,nodes=1,mem=3gb ${pathScripts}/bsub_script_cds
	
    fi
    
done

############################################################################################
