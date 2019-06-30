#!/bin/bash

export sp=$1
export type=$2
export i=$3

####################################################################################

export path=LncEvoDevo
export pathTools=Tools

export pathResults=${path}/results/mutual_information_network/${type}
export pathAracne=${pathTools}/ARACNe-AP/
export pathScripts=${path}/scripts/mutual_information_network

####################################################################################

if [ ${type} = "all_genes" ]; then
    export regulators=AllGenes.txt
fi

if [ ${type} = "lncRNAs_only" ]; then
    export regulators=LncRNAs.txt
fi

####################################################################################

echo "#!/bin/bash" >  ${pathScripts}/bsub_script_aracne

echo "java -Xmx8G -jar ${pathAracne}/dist/aracne.jar -e ${pathResults}/${sp}/TPM.txt -o ${pathResults}/${sp}/ --tfs ${pathResults}/${regulators}  --pvalue 1E-8 --threads 1 --seed ${i} --nodpi" >> ${pathScripts}/bsub_script_aracne

qsub -q huge -l s_rss=10G,sps=1 -o ${pathScripts}/std_output_aracne_${sp}.txt -e ${pathScripts}/std_error_aracne_${sp}.txt ${pathScripts}/bsub_script_aracne


####################################################################################
