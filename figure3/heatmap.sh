#!/bin/bash


dir=/analysis/lixf/CRISPRa_i/NCCIT/ChIP/H3K27ac/CRISPRa/

#################   N-DOX  L1PA   ################
FLL1="/LiuLab/reference/Human/GRCh38/TE/TE_bed/LINE-1.removechrY.FullLength.strand.bed"

############# cat enhancer.K562.hg38.bed|awk '{print $3-$2}'|datamash mean 1 median 1 q3 1
########## 2379.3844940745 1570    3040

########################  ########################
prefix=FLL1_combine_controlYAxis_removeK9me3
threads=70
mkdir -p ${dir}/4-Deeptools/H3K27ac/${prefix}/
cd ${dir}/4-Deeptools/H3K27ac/${prefix}/
# bedtools intersect -a ${FLL1} -b ${CpG} -F 0.5 -wa >CpG_FLL1.strand.bed
# bedtools intersect -a ${FLL1} -b ${CpG} -F 0.5 -wa -v >non_CpG_FLL1.strand.bed

# CpG_FL=CpG_FLL1.strand.bed
# non_CpG_FL=non_CpG_FLL1.strand.bed


########### H3K27ac ##################
{
computeMatrix scale-regions \
--regionBodyLength 6000 \
--startLabel "5'" \
--endLabel "3'" \
 -p ${threads} \
 --missingDataAsZero \
 -R ${FLL1} \
 --sortUsingSamples 7 \
 --binSize 50 \
 --skipZeros \
-S /analysis/lixf/CRISPRa_i/NCCIT/ChIP/H3K27ac/CRISPRa/4-bwCompare_bowtie2/random/H3K27ac_sgL1_vs_WT_subtract.bw \
/analysis/lixf/CRISPRa_i/NCCIT/ChIP/PolII_K9me3_RAD21_CTCF/4-bwCompare_bowtie2/random/PolII_A_sgL1_vs_WT_subtract.bw \
/analysis/lixf/CRISPRa_i/NCCIT/ChIP/H3K27ac/CRISPRi/4-bwCompare_bowtie2/random/H3K27ac_sgL1_vs_WT_subtract.bw \
/analysis/lixf/CRISPRa_i/NCCIT/ChIP/PolII_K9me3_RAD21_CTCF/4-bwCompare_bowtie2/random/PolII_I_sgL1_vs_WT_subtract.bw \
 -b 3000 -a 3000 \
 --samplesLabel K27ac_A PolII_A K27ac_I PolII_I \
 -o ${dir}/4-Deeptools/H3K27ac/${prefix}/H3K27ac.matrix.mat.gz  



plotHeatmap -m ${dir}/4-Deeptools/H3K27ac/${prefix}/H3K27ac.matrix.mat.gz \
            -o ${dir}/4-Deeptools/H3K27ac/${prefix}/H3K27ac_7.pdf \
            --outFileSortedRegions ${dir}/4-Deeptools/H3K27ac/${prefix}/H3K27ac.peak.bed \
           --heatmapHeight 15 --heatmapWidth 1.5 \
           --colorList "#2c4688,#C4CBDE,#FDEDED,#F9C0C2,#F05A61" "#2c4688,#C4CBDE,#FDEDED,#F9C0C2,#F05A61" "#2c4688,#556598,#7784AC,#C4CBDE,#FBD3D5" "#2c4688,#556598,#7784AC,#C4CBDE,#FBD3D5" \
            --colorNumber 300 \
            --regionsLabel  "FLL1"   \
            --startLabel "L" \
            --xAxisLabel " " \
            --endLabel "R"  \
            --whatToShow 'plot, heatmap and colorbar' \
            --sortUsingSamples 1 \
            --zMin -40 -30 -90 -90  --zMax 120 90 30 30  \
            --yMin -10 -10 -35 -45 --yMax 160 140 10 10 \
            --sortUsing mean   &



} &



wait




