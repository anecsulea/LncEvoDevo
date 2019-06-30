#!/bin/bash

export species=$1

#############################################################################

export path=LncEvoDevo

export pathDocs=${path}/docs
export pathAnnot=${path}/data/ensembl_annotations/${species}
export pathResults=${path}/results/stringtie_assembly/${species}
export pathScripts=${path}/scripts/transcript_assembly

export release=94

#############################################################################

if [ -e ${pathResults}/combined ]; then
    echo "path output already there"
else
    mkdir ${pathResults}/combined
fi

if [ -e ${pathResults}/combined/gtf_files.txt ]; then
    echo "gtf list already there, removing it"
    rm  ${pathResults}/combined/gtf_files.txt
fi

#############################################################################

for sample in `ls ${pathResults} | grep -v txt | grep -v combined`
do
    echo ${pathResults}/${sample}/filtered_transcripts.gtf >> ${pathResults}/combined/gtf_files.txt
done

#############################################################################

echo "#!/bin/bash" > bsub_script_stringtie

echo "stringtie -v --merge -G ${pathAnnot}/FilteredTranscripts_Ensembl${release}.gtf -m 150 -a 8 -p 8 -F 0 -T 0 -f 0.05 -o ${pathResults}/combined/assembled_transcripts.gtf ${pathResults}/combined/gtf_files.txt">> bsub_script_stringtie

qsub -q huge -l s_rss=10G,sps=1 -o ${pathScripts}/std_output_merge_stringtie_${species}.txt -e ${pathScripts}/std_error_merge_stringtie_${species}.txt ${pathScripts}/bsub_script_stringtie

#############################################################################
