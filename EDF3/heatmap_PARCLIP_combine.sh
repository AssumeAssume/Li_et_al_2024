#!/bin/bash

dir=/analysis2/lixf/DIS3/PAR_CLIP/


peaks=/LiuLab/reference/Human/GRCh38/TE/TE_bed/LINE-1.removechrY.FullLength.strand.bed

prefix=FLL1_mappablity_combine
threads=90
mkdir -p ${dir}/4-Deeptools/DIS3_OVERLAP/${prefix}/
cd ${dir}/4-Deeptools/DIS3_OVERLAP/${prefix}/


########### DIS3 ##################
{
computeMatrix scale-regions \
--regionBodyLength 6000 \
--startLabel "5'" \
--endLabel "3'" \
 -p ${threads} \
 --missingDataAsZero \
 -R ${peaks} \
 --sortUsingSamples 1 \
 --binSize 20 \
 --skipZeros \
-S /analysis2/lixf/DIS3/PAR_CLIP/4-bwCompare_bowtie2/random/PAR_CLIP_mean.bw \
/LiuLab/reference/Human/GRCh38/k100.Umap.MultiTrackMappability.bw \
-b 1000 -a 1000 \
 --samplesLabel PAR_DIS3  mappability \
 -o ${dir}/4-Deeptools/DIS3_OVERLAP/${prefix}/DIS3.matrix.mat.gz  

plotHeatmap -m ${dir}/4-Deeptools/DIS3_OVERLAP/${prefix}/DIS3.matrix.mat.gz \
            -o ${dir}/4-Deeptools/DIS3_OVERLAP/${prefix}/DIS3.pdf \
            --outFileSortedRegions ${dir}/4-Deeptools/DIS3_OVERLAP/${prefix}/DIS3.peak.bed \
           --heatmapHeight 15 --heatmapWidth 2.0 \
            --colorList "#2c4688,white,#f7aaac,#f16164,#eb1c25" \
            --colorNumber 300 \
            --regionsLabel "FL L1"  \
            --startLabel "L" \
            --xAxisLabel " " \
            --endLabel "R"  \
            --whatToShow 'heatmap and colorbar' \
            --sortUsingSamples 1 \
            --sortUsing mean   \
          --zMax 5 5 1 &

plotProfile -m ${dir}/4-Deeptools/DIS3_OVERLAP/${prefix}/DIS3.matrix.mat.gz \
--outFileNameData ${dir}/4-Deeptools/DIS3_OVERLAP/${prefix}/DIS3_profile_grouped.txt \
    -o profile.pdf


} &


wait




