#!/bin/sh


cd /analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/K562_WGS_INSERTION


fq=/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/K562_WGS_INSERTION/0-rawdata

ref=/LiuLab/reference/Human/GRCh38/LNC5/3TR_to_5TR_ont.mmi
# #gtfref=/LiuLab/reference/Human/GRCh38/GTF/gencode.v31.chr_patch_hapl_scaff.annotation.gff

# Contaminants=/LiuLab/reference/Contaminants/contaminant_list.txt
# ribokmers=/LiuLab/reference/Human/ribokmers.fa


mkdir 0-FastaQC

mkdir 2-Minimap2
mkdir 3-Sam2Bam



for i in `cat ./sh_file/file`; do
########## Quality Control ##########
# NanoStat -t 4 --fastq ${fq}/${i}.fq.gz -o ./0-FastaQC -p ${i}.nanostat -n ${i}.stat &
# NanoPlot -t 4 --fastq ${fq}/${i}.fq.gz --plots hex dot -o ./0-FastaQC -p ${i}.nanoplot &
############ mapping ################################################
minimap2 --sam-hit-only -t 72 -ax map-ont  $ref ${fq}/${i}.fq.gz > ./2-Minimap2/${i}.sam 
samtools view -bhSu -F 2304 -@ 72 ./2-Minimap2/${i}.sam > ./3-Sam2Bam/${i}.bam
samtools sort -@ 72 ./3-Sam2Bam/${i}.bam -o ./3-Sam2Bam/${i}_sorted.bam
samtools index -@ 72 ./3-Sam2Bam/${i}_sorted.bam ./3-Sam2Bam/${i}_sorted.bai 
samtools flagstat -@ 72 ./3-Sam2Bam/${i}_sorted.bam > ./0-FastaQC/${i}_sorted.stat
############ bigwig ################################################


done
