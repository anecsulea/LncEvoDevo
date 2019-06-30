#!/bin/bash

export sp=$1
export annot=$2

################################################################

export path=LncEvoDevo

export pathDocs=${path}/docs
export pathHisat=${path}/results/hisat/${sp}
export pathResults=${path}/results/expression_estimation/${sp}
export pathScripts=${path}/scripts/expression_estimation

################################################################

export release=94

################################################################

echo "#!/bin/bash" >  ${pathScripts}/qsub_script

export running="0"

for sample in `ls ${pathHisat} | grep -v txt`
do

    if [ -e ${pathResults}/${sample}/ReadCounts_${annot}.txt ]; then
	echo ${sample} "already done"
    else
	export running="1"
	
	if [ ${sp} = "Chicken" ]; then
	    export library=`grep ^${sample}$'\t' ${pathDocs}/LibraryType_Chicken.txt | cut -f 2`
	else
	    export library=fr-firststrand
	fi
	
	export strand="2"
	
	if [ ${library} = "fr-unstranded" ]; then
	    export strand="0"
	fi
	
	if [ ${library} = "fr-secondstrand" ]; then
	    export strand="1"
	fi
	
	echo ${sample} ${library} ${strand}
	
	echo "Rscript --vanilla ${pathScripts}/count.reads.subread.R in2p3 ${annot} ${release} ${sp} ${sample} 1 ${strand}" >> ${pathScripts}/qsub_script
    fi    
done

if [ ${running} = "1" ]; then
  
    qsub -q long -l s_rss=4G,sps=1 -o ${pathScripts}/std_output_subread_${sp}_${annot}.txt -e ${pathScripts}/std_error_subread_${sp}_${annot}.txt ${pathScripts}/qsub_script
   
fi

################################################################
