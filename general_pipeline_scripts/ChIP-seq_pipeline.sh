#!/bin/bash

bowtie2ref="/LiuLab/reference/Human/GRCh38/bowtie2/GRCh38.p12.genome"


data_grep="$"

metadata="/analysis/lixf/CTBP1/ChIP/human/K562//PolII_H3K27ac.metadata.txt"

outdir="/analysis/lixf/CTBP1/ChIP/human/K562/"


thread=6

cd ${outdir}

# get srr list and metadata information
cat ${metadata}|grep -Pi "${data_grep}"|awk '{print $2}' >SRR.list
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$2]=$0}NR>FNR{print a[$1]}' ${metadata}  SRR.list > metadata.txt

declare -A array

eval $(cat ${metadata} |grep -Pi "${data_grep}"|sed '/^$/d'|awk -F "\t" '{print "array["$1"]="$2}')


mkdir -p ${outdir}/0-QC/flag_stat
mkdir -p ${outdir}/0-QC/clean_qc 
mkdir -p ${outdir}/0-clean_data/paired
mkdir -p ${outdir}/0-clean_data/unpaired
mkdir -p ${outdir}/0-clean_data/no_ribo
mkdir -p ${outdir}/0-clean_data/single_end

mkdir -p ${outdir}/1-bowtie2/sam
mkdir -p ${outdir}/1-bowtie2/bam/uniq
mkdir -p ${outdir}/1-bowtie2/bam/random
mkdir -p ${outdir}/1-bowtie2/bw/uniq
mkdir -p ${outdir}/1-bowtie2/bw/random
mkdir -p ${outdir}/3-MACS2_bowtie2/
mkdir -p ${outdir}/4-Deeptools/


mkdir -p ${outdir}/7-TagDirectory_bowtie2/uniq
mkdir -p ${outdir}/8-analyzeRepeats_bowtie2/uniq

mkdir -p ${outdir}/7-TagDirectory_bowtie2/random
mkdir -p ${outdir}/8-analyzeRepeats_bowtie2/random
# clean data

for key in ${!array[@]}
    do
    
    {

        fastqc ${array[$key]}* -t 4 -o ${outdir}/0-QC &

        trim_galore -j ${thread} -q 20 --phred33 \
            --no_report_file --trim-n --paired \
            ${array[$key]}_*1.fq.gz  ${array[$key]}_*2.fq.gz \
            -o ${outdir}//0-clean_data/paired/ \
            --basename ${key}

                echo "${key} cleaned"

        fastqc ${outdir}//0-clean_data/paired/${key}* -t 4 -o ${outdir}/0-QC/  &


        bowtie2 -p ${thread} --no-mixed --no-discordant \
            --end-to-end --very-sensitive --maxins 700  \
            -x ${bowtie2ref} \
            -1 ${outdir}//0-clean_data/paired/${key}_val_1.fq.gz \
            -2 ${outdir}//0-clean_data/paired/${key}_val_2.fq.gz | samtools view -bhSu  | samtools sort -@ ${thread} -m 4G -O bam -T ${key}_q1.prefix |samtools rmdup - ${outdir}/1-bowtie2/bam/random/${key}_random.sorted.rudup.bam

        samtools index ${outdir}/1-bowtie2/bam/random/${key}_random.sorted.rudup.bam
        samtools flagstat ${outdir}/1-bowtie2/bam/random/${key}_random.sorted.rudup.bam > ${outdir}//0-QC/flag_stat/${key}_random.sorted.stat

        bamCoverage -b ${outdir}/1-bowtie2/bam/random/${key}_random.sorted.rudup.bam \
                    --normalizeUsing RPKM \
                    -p ${thread} \
                    --binSize 5 \
                    -o ${outdir}/1-bowtie2/bw/random/${key}.normalized.bw \
                    --ignoreDuplicates \
                    --extendReads
        ##### extract uniq mapping reads #####
        #### filter MAPQ >10 READS ####
        samtools view -b ${outdir}/1-bowtie2/bam/random/${key}_random.sorted.rudup.bam -q 10 > ${outdir}/1-bowtie2/bam/uniq/${key}_uniq.sorted.rudup.bam

        samtools index ${outdir}/1-bowtie2/bam/uniq/${key}_uniq.sorted.rudup.bam
        samtools flagstat ${outdir}/1-bowtie2/bam/uniq/${key}_uniq.sorted.rudup.bam > ${outdir}//0-QC/flag_stat/${key}_uniq.sorted.stat

        bamCoverage -b ${outdir}/1-bowtie2/bam/uniq/${key}_uniq.sorted.rudup.bam \
                    --normalizeUsing RPKM \
                    -p ${thread} \
                    --binSize 5 \
                    -o ${outdir}/1-bowtie2/bw/uniq/${key}.normalized.bw \
                    --ignoreDuplicates \
                    --extendReads
###### random

        makeTagDirectory ${outdir}/7-TagDirectory_bowtie2/random/${key}  ${outdir}/1-bowtie2/bam/random/${key}_random.sorted.rudup.bam

        analyzeRepeats.pl repeats hg38 -norm 1e7 -d ${outdir}/7-TagDirectory_bowtie2/random/${key} -noCondensing -strand both -upstream 1000 -downstream 1000 >${outdir}/8-analyzeRepeats_bowtie2/random/${key}.norm.noCondensing.txt 
        analyzeRepeats.pl repeats hg38 -raw -d ${outdir}/7-TagDirectory_bowtie2/random/${key} -noCondensing -strand both -upstream 1000 -downstream 1000 >${outdir}/8-analyzeRepeats_bowtie2/random/${key}.raw.noCondensing.txt 

        analyzeRepeats.pl repeats hg38 -norm 1e7 -d ${outdir}/7-TagDirectory_bowtie2/random/${key} -condenseL1 -strand both -upstream 1000 -downstream 1000 >${outdir}/8-analyzeRepeats_bowtie2/random/${key}.norm.subfamily.txt 
        analyzeRepeats.pl repeats hg38 -raw -d ${outdir}/7-TagDirectory_bowtie2/random/${key} -condenseL1 -strand both -upstream 1000 -downstream 1000 >${outdir}/8-analyzeRepeats_bowtie2/random/${key}.raw.subfamily.txt 
###### uniq 

        makeTagDirectory ${outdir}/7-TagDirectory_bowtie2/uniq/${key}  ${outdir}/1-bowtie2/bam/uniq/${key}_uniq.sorted.rudup.bam

        analyzeRepeats.pl repeats hg38 -norm 1e7 -d ${outdir}/7-TagDirectory_bowtie2/uniq/${key} -noCondensing -strand both -upstream 1000 -downstream 1000 >${outdir}/8-analyzeRepeats_bowtie2/uniq/${key}.norm.noCondensing.txt 
        analyzeRepeats.pl repeats hg38 -raw -d ${outdir}/7-TagDirectory_bowtie2/uniq/${key} -noCondensing -strand both -upstream 1000 -downstream 1000 >${outdir}/8-analyzeRepeats_bowtie2/uniq/${key}.raw.noCondensing.txt 

        analyzeRepeats.pl repeats hg38 -norm 1e7 -d ${outdir}/7-TagDirectory_bowtie2/uniq/${key} -condenseL1 -strand both -upstream 1000 -downstream 1000 >${outdir}/8-analyzeRepeats_bowtie2/uniq/${key}.norm.subfamily.txt 
        analyzeRepeats.pl repeats hg38 -raw -d ${outdir}/7-TagDirectory_bowtie2/uniq/${key} -condenseL1 -strand both -upstream 1000 -downstream 1000 >${outdir}/8-analyzeRepeats_bowtie2/uniq/${key}.raw.subfamily.txt 


    }  &

done

wait

cd ${outdir}/
multiqc . -o ${outdir}//0-QC


