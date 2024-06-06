#!/bin/bash

dir=/analysis/lixf/sgRNA/human/LINE_screening_combine/

#################     ################
FLL1=/LiuLab/reference/Human/GRCh38/TE/TE_bed/LINE-1.removechrY.FullLength.strand.bed
CpG=/LiuLab/reference/Human/GRCh38/CpG_island/unmasked_CpG.bed

FLL1PA1=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1HS.bed
FLL1PA2=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1PA2.bed
FLL1PA3=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1PA3.bed
FLL1PA4=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1PA4.bed
FLL1PA5=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1PA5.bed
FLL1PA6=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1PA6.bed
FLL1PA7=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1PA7.bed
FL_others=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_NonL1PA1_7_LINE1.strand.bed

########################  ########################
prefix=FLL1PA1234567_WT
threads=70
mkdir -p ${dir}/4-Deeptools/FLL1_profile/${prefix}/
cd ${dir}/4-Deeptools/FLL1_profile/${prefix}/



###########  ##################
{
computeMatrix scale-regions \
 --regionBodyLength 6000 \
 -p ${threads} \
 --startLabel "5'" \
 --endLabel "3'" \
 --missingDataAsZero \
 -R ${FLL1PA1} ${FLL1PA2} ${FLL1PA3} ${FLL1PA4} ${FLL1PA5} ${FLL1PA6} ${FLL1PA7} ${FL_others} \
 --sortUsingSamples 1 \
 --binSize 100 \
 --skipZeros \
-S  /analysis/lixf/CTBP1/ChIP/human/K562//4-bwCompare_bowtie2/random/H3K27ac_WT_mean.bw \
/analysis/lixf/CTBP1/ChIP/human/K562//4-bwCompare_bowtie2/random/H3K9me3_WT_mean.bw \
 -b 3000 -a 3000 \
 --samplesLabel  K27ac K9me3 \
 -o ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine.matrix.mat.gz  



# plotHeatmap -m ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine.matrix.mat.gz \
#             -o ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine_bilinear_1.pdf \
#             --outFileSortedRegions ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine.peak.bed \
#            --heatmapHeight 15 --heatmapWidth 2.5 \
#            --interpolationMethod "bilinear" \
#            --colorList "#2c4688,#7889B3,white,#eb1c25"  "#2c4688,#7889B3,white,#eb1c25" "#2c4688,white,#f7aaac,#eb1c25" "#2c4688,white,#eb1c25"  \
#             --colorNumber 300 \
#             --startLabel "5'" \
#             --xAxisLabel " " \
#             --endLabel "3'"  \
#             --whatToShow 'heatmap and colorbar' \
#             --sortUsingSamples 2 \
#             --regionsLabel "PA1" "PA2" "PA3" "PA4" "PA5" "PA6" "PA7" "others" \
#             --zMax 100 3 3 50 \
#             --sortUsing mean   &


plotHeatmap -m ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine.matrix.mat.gz \
            -o ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine_bilinear_2.pdf \
            --outFileSortedRegions ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine.peak.bed \
           --heatmapHeight 20 --heatmapWidth 3.0 \
           --interpolationMethod "bilinear" \
           --colorList "#2c4688,white,#f7aaac,#eb1c25"  "#2c4688,white,#eb1c25" "#2c4688,white,#eb1c25" "#2c4688,white,#eb1c25" "#2c4688,white,#eb1c25" "#2c4688,white,#eb1c25" "#2c4688,white,#f7aaac,#eb1c25" "#2c4688,white,#f7aaac,#eb1c25" "#2c4688,white,#f7aaac,#eb1c25" "#2c4688,white,#f7aaac,#eb1c25" "#2c4688,white,#f7aaac,#eb1c25"  "#2c4688,white,#f7aaac,#eb1c25" "#2c4688,white,#eb1c25"  \
            --colorNumber 300 \
            --startLabel "5'" \
            --xAxisLabel " " \
            --endLabel "3'"  \
            --whatToShow 'heatmap and colorbar' \
            --sortUsingSamples 1 \
            --regionsLabel "PA1" "PA2" "PA3" "PA4" "PA5" "PA6" "PA7" "others" \
            --zMax 2 100 100 100 100 100 2 3 2.5 2 5 20 50 \
            --sortUsing mean   &




plotProfile -m ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine.matrix.mat.gz \
              -out ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine_profile_grouped.pdf \
              --outFileNameData ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine_profile_grouped.txt \
              --numPlotsPerRow 8 \
              --plotWidth 5 \
              --plotHeight 6 \
              --yMax 50 130 2 330 3 100 140 \
            --yMin 0 0 0 0 0 0 0  \
              --legendLocation upper-right 
             # --colors "#bd3106" "#d9700e" "#eebe04" "#c3d6ce" "#89a6bb" "#454b87" & 



### --colorList "#2c4688,white,#f7aaac,#f16164,#eb1c25" \
} &

wait




