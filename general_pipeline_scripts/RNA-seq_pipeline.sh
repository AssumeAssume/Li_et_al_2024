#!/bin/bash

ref="/LiuLab/reference/Human/GRCh38/hisat2/GRCh38.p12.genome.fa"


data_grep="CTBP1"

metadata="/analysis/lixf/RNA/human/K562/KO/metadata.txt"

outdir="/analysis/lixf/RNA/human/K562/KO/CTBP1/"


KnownGeneref=/LiuLab/reference/Human/GRCh38/GTF/gencode.v31.chr_patch_hapl_scaff.annotation.gff
RepeatMasker=/LiuLab/reference/Human/GRCh38/TE/Homo_sapiens.hg38.UCSC.RepeatMasker.gtf
RepeatMaskerUnique=/LiuLab/reference/Human/GRCh38/TE/Homo_sapiens.hg38.UCSC.RepeatMasker.unique.gtf

thread=24

cd ${outdir}

# get srr list and metadata information
cat ${metadata}|grep -Pi "${data_grep}"|awk '{print $2}' >SRR.list
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$2]=$0}NR>FNR{print a[$1]}' ${metadata}  SRR.list >metadata.txt

declare -A array

eval $(cat ${metadata} |grep -Pi "${data_grep}"|sed '/^$/d'|awk -F "\t" '{print "array["$1"]="$2}')

mkdir -p ${outdir}/0-QC/flag_stat
mkdir -p ${outdir}/0-clean_data/paired
mkdir -p ${outdir}/0-clean_data/unpaired
mkdir -p ${outdir}/0-clean_data/no_ribo
mkdir -p ${outdir}/0-clean_data/single_end
mkdir -p ${outdir}/1-fastq2sam2bam/bam
mkdir -p ${outdir}/1-fastq2sam2bam/sam
mkdir -p ${outdir}/2-bam2BigWig/
mkdir -p ${outdir}/3-featureCounts/KnownGene
mkdir -p ${outdir}/3-featureCounts/RepeatMasker
mkdir -p ${outdir}/3-featureCounts/RepeatMasker_Unique
mkdir -p ${outdir}/4-RawCount/KnownGene
mkdir -p ${outdir}/4-RawCount/RepeatMasker
mkdir -p ${outdir}/4-RawCount/RepeatMasker_Unique
mkdir -p ${outdir}/5-DESeq2/ALL
mkdir -p ${outdir}/5-DESeq2/KnownGene
mkdir -p ${outdir}/5-DESeq2/RepeatMasker
mkdir -p ${outdir}/5-DESeq2/RepeatMasker_Unique



# clean data

for key in ${!array[@]}
    do

        {
fastqc ${array[$key]}_* -t 4 -o ${outdir}/0-QC &

            trim_galore -j ${thread} -q 20 --phred33 \
            --no_report_file --trim-n --paired \
            ${array[$key]}_*R1*fastq.gz  ${array[$key]}_*R2*fastq.gz \
            -o ${outdir}//0-clean_data/paired/ \
            --basename ${key}

            echo "${key} cleaned"

fastqc ${outdir}//0-clean_data/paired/${key}* -t 4 -o ${outdir}/0-QC/  &
# mapping
hisat2 -p ${thread} --dta --no-mixed --no-discordant \
    -x $ref \
    --rna-strandness RF \
    -1 ${outdir}//0-clean_data/paired/${key}_val_1.fq.gz \
    -2  ${outdir}//0-clean_data/paired/${key}_val_2.fq.gz \
    -S ${outdir}//1-fastq2sam2bam/sam/${key}.sam
    samtools view \
         -bhSu ${outdir}//1-fastq2sam2bam/sam/${key}.sam | \
        samtools sort -@ ${thread} -m 4G -O bam -T out.prefix  -o ${outdir}//1-fastq2sam2bam/bam/${key}.sorted.bam
samtools index ${outdir}//1-fastq2sam2bam/bam/${key}.sorted.bam
samtools flagstat ${outdir}//1-fastq2sam2bam/bam/${key}.sorted.bam > ${outdir}//0-QC/flag_stat/${key}.sorted.stat
# bam2bigwig
SIZE=$(awk 'NR==5{print $1}' ${outdir}//0-QC/flag_stat/${key}.sorted.stat)
SCALE=$(echo "scale=8;1000000/$SIZE" | bc)
bamCoverage -b ${outdir}//1-fastq2sam2bam/bam/${key}.sorted.bam \
            --scaleFactor $SCALE \
            -p ${thread} \
            --binSize 10 \
            -o ${outdir}//2-bam2BigWig/${key}.normalized.bw


# featurecounts

featureCounts -T ${thread} -p -t exon -g gene_name -Q 10 -s 2 -a ${KnownGeneref} -o ${outdir}/3-featureCounts/KnownGene/${key}.KnownGene.counts.txt ${outdir}/1-fastq2sam2bam/bam/${key}.sorted.bam
awk '{if($0 !~ "^#"){print $1"\t"$7}}' ${outdir}/3-featureCounts/KnownGene/${key}.KnownGene.counts.txt > ${outdir}/4-RawCount/KnownGene/${key}.KnownGene.rawcounts.tsv 

featureCounts -T ${thread} -p -t exon -g gene_id -Q 0 -O -M --fraction \
               -a $RepeatMasker -o ${outdir}/3-featureCounts/RepeatMasker/${key}.RepeatMasker.counts.txt ${outdir}/1-fastq2sam2bam/bam/${key}.sorted.bam
awk '{if($0 !~ "^#"){print $1"\t"$7}}' ${outdir}/3-featureCounts/RepeatMasker/${key}.RepeatMasker.counts.txt > ${outdir}/4-RawCount/RepeatMasker/${key}.RepeatMasker.rawcounts.tsv 


featureCounts -T ${thread} -p -t exon -g unique_id -f -Q 10 \
               -a $RepeatMaskerUnique -o ${outdir}/3-featureCounts/RepeatMasker_Unique/${key}.RepeatMasker.unique.counts.txt ${outdir}/1-fastq2sam2bam/bam/${key}.sorted.bam
awk '{if($0 !~ "^#"){print $1"\t"$7}}' ${outdir}/3-featureCounts/RepeatMasker_Unique/${key}.RepeatMasker.unique.counts.txt > ${outdir}/4-RawCount/RepeatMasker_Unique/${key}.RepeatMasker.unique.rawcounts.tsv 




             } & 



done
wait



cd ${outdir}//0-QC
multiqc .

