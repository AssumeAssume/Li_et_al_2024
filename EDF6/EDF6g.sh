#!/bin/bash

dir=/analysis/lixf/sgRNA/human/LINE_screening_combine/figure2/KMT2D


FLL1=/LiuLab/reference/Human/GRCh38/TE/TE_bed/LINE-1.removechrY.FullLength.strand.bed
CpG=/LiuLab/reference/Human/GRCh38/CpG_island/unmasked_CpG.bed

DEL1=/analysis/lixf/RNA/human/K562/KO/KMT2D/6-LINE1M/intron_inter_DEL1.bed


##
prefix=FLL1_separate_K4me1_K27ac
threads=70
mkdir -p ${dir}/4-Deeptools/FLL1_profile/${prefix}/
cd ${dir}/4-Deeptools/FLL1_profile/${prefix}/



{
computeMatrix scale-regions \
 --regionBodyLength 6000 \
 -p ${threads} \
 --startLabel "5'" \
 --endLabel "3'" \
 --missingDataAsZero \
 -R ${FLL1} \
 --sortUsingSamples 1 \
 --binSize 100 \
 --skipZeros \
-S  /analysis/lixf/KMT2D/ChIP/human/K562/4-bwCompare_bowtie2/H3K4me1_WT_mean.bw \
/analysis/lixf/KMT2D/ChIP/human/K562/4-bwCompare_bowtie2/H3K4me1_KO_mean.bw \
/analysis/lixf/KMT2D/ChIP/human/K562/WY_20210913_H3K27ac/4-bwCompare_bowtie2/H3K27ac_WT_mean.bw \
 /analysis/lixf/KMT2D/ChIP/human/K562/WY_20210913_H3K27ac/4-bwCompare_bowtie2/H3K27ac_KO_mean.bw \
 -b 6000 -a 6000 \
 --samplesLabel K4me1_WT K4me1_KO K27ac_WT  K27ac_KO \
 -o ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine.matrix.mat.gz  




plotHeatmap -m ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine.matrix.mat.gz \
            -o ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine_bilinear_2.pdf \
            --outFileSortedRegions ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine.peak.bed \
           --heatmapHeight 20 --heatmapWidth 3 \
           --interpolationMethod "bilinear" \
            --colorList "#2c4688,#475E97,#5F73A5,white,#F48A90,#F38187,#F1686F,#EE4951,#eb1c25" "#2c4688,#475E97,#5F73A5,white,#F48A90,#F38187,#F1686F,#EE4951,#eb1c25" \
            --colorNumber 300 \
            --regionsLabel "FLL1"  \
            --startLabel "5'" \
            --xAxisLabel " " \
            --endLabel "3'"  \
            --whatToShow 'heatmap and colorbar' \
            --sortUsingSamples 1 \
            --zMin 0 0 0 0 --zMax 50 50 20 20 \
            --sortUsing mean   &

} &

wait




