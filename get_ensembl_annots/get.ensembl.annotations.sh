#######################################################################

export species=$1
export release=$2

#######################################################################

if [ ${species} = "Mouse" ]; then
    export db="mus_musculus_core_${release}_38"
fi

if [ ${species} = "Rat" ]; then
    if [ "$release" -ge 80 ]; then
	echo "using rn6"
	export db="rattus_norvegicus_core_${release}_6"
    else
	echo "using rn5"
	export db="rattus_norvegicus_core_${release}_5"
    fi
fi

if [ ${species} = "Human" ]; then
    export db="homo_sapiens_core_${release}_38"
fi


if [ ${species} = "Macaque" ]; then
    export db="macaca_mulatta_core_${release}_10"
fi


if [ ${species} = "Chicken" ]; then
    if [ "$release" -ge 86 ]; then
	export db="gallus_gallus_core_${release}_5"
    else
	export db="gallus_gallus_core_${release}_4"
    fi
fi


#######################################################################

if [ -e ../../data/ensembl_annotations/${species} ]; then 
    echo "path exists"
else
    mkdir ../../data/ensembl_annotations/${species}
fi

echo "extracting ensembl_annotations for "${species}

#######################################################################

## exon coordinates

echo "use $db; " > get.exons.${species}.sh
echo "select exon.stable_id, seq_region.name, exon.seq_region_start, exon.seq_region_end, exon.seq_region_strand from exon, seq_region where exon.seq_region_id=seq_region.seq_region_id; " >> get.exons.${species}.sh

mysql -h ensembldb.ensembl.org -P 5306 -u anonymous < get.exons.${species}.sh > ../../data/ensembl_annotations/${species}/ExonCoords_Ensembl${release}.txt

## transcript info

echo "use $db; " > get.transcripts.${species}.sh
echo "select gene.stable_id, transcript.stable_id, transcript.biotype, seq_region.name, transcript.seq_region_start, transcript.seq_region_end, transcript.seq_region_strand from transcript, gene, seq_region where transcript.gene_id=gene.gene_id and transcript.seq_region_id=seq_region.seq_region_id; " >> get.transcripts.${species}.sh

mysql -h ensembldb.ensembl.org -P 5306 -u anonymous < get.transcripts.${species}.sh > ../../data/ensembl_annotations/${species}/TranscriptInfo_Ensembl${release}.txt

## gene info 

echo "use $db; " > get.genes.${species}.sh
echo "select gene.stable_id, gene.biotype, gene.status, gene.description, seq_region.name, gene.seq_region_start, gene.seq_region_end, gene.seq_region_strand from gene, seq_region where gene.seq_region_id=seq_region.seq_region_id; " >> get.genes.${species}.sh

mysql -h ensembldb.ensembl.org -P 5306 -u anonymous < get.genes.${species}.sh > ../../data/ensembl_annotations/${species}/GeneInfo_Ensembl${release}.txt

## exons-transcripts

echo "use $db; " > get.exontx.${species}.sh
echo "select exon.stable_id, transcript.stable_id from exon_transcript, exon, transcript where exon_transcript.exon_id=exon.exon_id and exon_transcript.transcript_id=transcript.transcript_id; " >> get.exontx.${species}.sh

mysql -h ensembldb.ensembl.org -P 5306 -u anonymous < get.exontx.${species}.sh > ../../data/ensembl_annotations/${species}/ExonsTranscripts_Ensembl${release}.txt


mv *${species}.sh species_scripts/

#######################################################################
