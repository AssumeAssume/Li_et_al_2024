#!/bin/bash

###cd /analysis2/hongyq/Hi-C/CRISPRa_i/run1_2


Juicer=/analysis2/hongyq/script/juicer.sh
ref=/analysis2/reference/Human/GRCh38_primary/GRCh38.primary_assembly.genome.fa
site="MboI"
chromSizes=/analysis2/reference/Human/GRCh38_primary/GRCh38.primary_assembly.genome.chrom.sizes
restrictionSite=/analysis2/reference/Human/GRCh38_primary/restriction_sites/GRCh38.primary_assembly.genome_${site}.txt
scriptDir=/analysis2/software/juicer


threads=50
genomeID=GRCh38


metadata=/analysis2/lixf/Hi-C/human/K562/BLY_202320214_CTBP1/hic_metadata.txt
outdir="/analysis2/lixf/Hi-C/human/K562/BLY_202320214_CTBP1/"
cd ${outdir}

# get srr list and metadata information
cat ${metadata}|grep -Pi "${data_grep}"|awk '{print $2}' >SRR.list
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$2]=$0}NR>FNR{print a[$1]}' ${metadata}  SRR.list > metadata.txt

declare -A array

eval $(cat ${metadata} |grep -Pi "${data_grep}"|sed '/^$/d'|datamash -s -g 1 collapse 2 |awk -F "\t" '{print "array["$1"]="$2}')

mkdir -p  ${outdir}/1-Juicer/


for key in ${!array[@]}
	do
	{
		topDir=${outdir}/1-Juicer/$key/
		
		mkdir -p ${topDir}/fastq
		readsarray=(`echo ${array[$key]}|tr ',' ' '`)
		for read in ${readsarray[@]}
			do
			ln -s ${read}_*1.fq.gz  ${topDir}/fastq/${read##*/}_R1.fastq.gz
			ln -s ${read}_*2.fq.gz  ${topDir}/fastq/${read##*/}_R2.fastq.gz
		done

		$Juicer -z $ref -d $topDir -s $site -S "early" -p $chromSizes -y $restrictionSite -D $scriptDir -t $threads -g $genomeID
		$Juicer -z $ref -d $topDir -s $site -S "final" -p $chromSizes -y $restrictionSite -D $scriptDir -t $threads -g $genomeID
# $Juicer -z $ref -d $topDir -s $site -S "postproc" -p $chromSizes -y $restrictionSite -D $scriptDir -t $threads -g $genomeID
	} &
done