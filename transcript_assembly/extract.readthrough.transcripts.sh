#!/bin/bash

export species=$1
export type=$2

##############################################################

export path=LncEvoDevo

export pathAnnot=${path}/data/ensembl_annotations/${species}
export pathResults=${path}/results/${type}_assembly/${species}
export pathScripts=${path}/scripts/transcript_assembly

if [ ${type} = "cufflinks" ]; then
    export filename="transcripts"
else
    if [ ${type} = "stringtie" ]; then
	export filename="assembled_transcripts"
    else
	echo "unknown type"
	exit
    fi
fi

export release=94

###########################################################################

for sample in `ls ${pathResults} | grep -v txt `
do

    if [ -e ${pathResults}/${sample}/readthrough_transcripts.txt ]; then
	echo "already done"
    else

	echo "#!/bin/bash" >  ${pathScripts}/bsub_script_readthrough

	echo "perl ${pathScripts}/extract.readthrough.transcripts.pl --pathEnsemblGTF=${pathAnnot}/FilteredTranscripts_Ensembl${release}.gtf --pathAssembledGTF=${pathResults}/${sample}/${filename}.gtf --pathGeneInfo=${pathAnnot}/GeneInfo_Ensembl${release}.txt --monoexonicBiotypes=\"Mt_tRNA,Mt_rRNA,IG_D_gene,IG_J_gene,snoRNA,misc_RNA,miRNA,snRNA,rRNA\" --pathOutput=${pathResults}/${sample}/readthrough_transcripts.txt"  >>  ${pathScripts}/bsub_script_readthrough
	

	qsub -q long -l s_rss=4G,sps=1 -o ${pathScripts}/std_output_readthrough_${species}.txt -e ${pathScripts}/std_error_readthrough_${species}.txt ${pathScripts}/bsub_script_readthrough
	
    fi
done

###########################################################################

