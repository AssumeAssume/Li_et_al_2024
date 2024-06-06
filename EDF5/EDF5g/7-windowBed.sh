#!/bin/bash

dir=$1

###### ls ${dir}/3-MACS2|tr '\n' ' '
prefix=$2
bed=$3
ymin=$4
ymax=$5
mkdir -p ${dir}/7-windowBed/${prefix}


queryBed=/analysis/lixf/sgRNA/human/LINE_screening_combine/CTBP1_KMT2D_fusion_transcript_analysis/2-taco_output/anno/TSS_anno.nameCol4.sorted.bed

sizeList=(0 1 2 4 8 16 32 64 128 256 512 1024 2048 )
echo $bed
echo $queryBed
outdir=${dir}/7-windowBed/${prefix}
for size in ${sizeList[@]}; do
{
	size_kb=`echo "${size}*1000"|bc -l`
	echo ${size_kb}

	bedtools window -a $bed -b $queryBed -l ${size_kb} -r 0 -sw >${outdir}/overlap_${size}K_upstream.bed
	bedtools window -a $bed -b $queryBed -r ${size_kb} -l 0 -sw >${outdir}/overlap_${size}K_downstream.bed

	cat ${outdir}/overlap_${size}K_upstream.bed |cut -f 10 |sort|uniq |awk -v distance="$size"  '{print $1,-distance}' > ${outdir}/overlap_${size}K_upstream_genename.txt

	cat ${outdir}/overlap_${size}K_downstream.bed |cut -f 10 |sort|uniq |awk -v distance="$size"  '{print $1,distance}' > ${outdir}/overlap_${size}K_downstream_genename.txt

} &
done

#sizeList2=(0 20 40 60 80 100 120 140 160 180 200)
#for size in ${sizeList2[@]}; do
#{
#	size_kb=`echo "${size}*1000"|bc -l`
#	echo ${size_kb}
#
##	bedtools window -a $bed -b $queryBed -l ${size_kb} -r 0 -sw >${outdir}/overlap_${size}K_upstream_arithmetic.bed
#	bedtools window -a $bed -b $queryBed -r ${size_kb} -l 0 -sw >${outdir}/overlap_${size}K_downstream_arithmetic.bed
#
#	cat ${outdir}/overlap_${size}K_upstream_arithmetic.bed |cut -f 10 |sort|uniq |awk -v distance="$size"  '{print $1,-distance}' > ${outdir}/overlap_${size}K_upstream_arithmetic_genename.txt

#	cat ${outdir}/overlap_${size}K_downstream_arithmetic.bed |cut -f 10 |sort|uniq |awk -v distance="$size"  '{print $1,distance}' > ${outdir}/overlap_${size}K_downstream_arithmetic_genename.txt

#} &
#done
wait
#/analysis/lixf/sgRNA/human/LINE_screening_combine/scripts/figure4_CTBP1/7-nerar_expression_plot.R $dir $prefix $ymin $ymax

