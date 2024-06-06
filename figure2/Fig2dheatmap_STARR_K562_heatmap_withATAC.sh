#!/bin/bash

dir=/analysis/lixf/sgRNA/human/LINE_screening_combine/STARR_related

#################   N-DOX  L1PA   ################
FLL1=/LiuLab/reference/Human/GRCh38/TE/TE_bed/LINE-1.removechrY.FullLength.strand.bed
CpG=/LiuLab/reference/Human/GRCh38/CpG_island/unmasked_CpG.bed

FLL1M=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1_separate/FLL1M.bed
FLL1PA1=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1_separate/FLL1PA1.bed
FLL1PA2=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1_separate/FLL1PA2.bed
FLL1PA3=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1_separate/FLL1PA3.bed
FLL1PA4=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1_separate/FLL1PA4.bed
FLL1PA5=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1_separate/FLL1PA5.bed
FLL1PA6=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1_separate/FLL1PA6.bed
FLL1PA7=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1_separate/FLL1PA7.bed
FLL1PA8_17=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1_separate/FLL1PA8_17.bed
FLL1PB=/LiuLab/reference/Human/GRCh38/TE/TE_bed/FL_L1_separate/FLL1PB.bed

######################## N_LF_KDM2B ########################
prefix=FLL1PA1234567_WT_with_ATAC
threads=70
mkdir -p ${dir}/4-Deeptools/FLL1_profile/${prefix}/
cd ${dir}/4-Deeptools/FLL1_profile/${prefix}/



########### H3K4me3 ##################
{
computeMatrix scale-regions \
 --regionBodyLength 6000 \
 -p ${threads} \
 --startLabel "5'" \
 --endLabel "3'" \
 --missingDataAsZero \
 -R ${FLL1PA1} ${FLL1PA2} ${FLL1PA3} ${FLL1PA4} ${FLL1PA5} ${FLL1PA6} ${FLL1PA7} ${FLL1PA8_17} ${FLL1PB} ${FLL1M}  \
 --sortUsingSamples 1 \
 --binSize 100 \
 --skipZeros \
-S  /analysis/lixf/tracks/STARR/several_cellLine/rawdata_analysis/StarrPeaker/K562_STARR.fc.bw \
/analysis/lixf/tracks/STARR/several_cellLine/rawdata_analysis/StarrPeaker/SHSY5Y_STARR.fc.bw \
/analysis/lixf/tracks/STARR/several_cellLine/rawdata_analysis/StarrPeaker/HCT116_STARR.fc.bw \
/analysis/lixf/tracks/STARR/several_cellLine/rawdata_analysis/StarrPeaker/HepG2_STARR.fc.bw \
/analysis/lixf/tracks/STARR/several_cellLine/rawdata_analysis/StarrPeaker/MCF7_STARR.fc.bw \
/analysis/lixf/tracks/STARR/several_cellLine/rawdata_analysis/StarrPeaker/A549_STARR.fc.bw  \
/analysis/lixf/CTBP1/ChIP/human/K562/WY_20210625_PolII_H3K27ac/4-bwCompare_bowtie2/random/H3K27ac_WT_mean.bw \
/analysis/lixf/CTBP1/ChIP/human/K562/WY_20220328_H3K4me1/4-bwCompare_bowtie2/random/H3K4me1_WT_mean.bw \
/analysis/lixf/CTBP1/ChIP/human/K562/WY_20210831_H3K4me3_K9Ac/4-bwCompare_bowtie2/random/H3K4me3_WT_mean.bw \
/analysis/lixf/CTBP1/ChIP/human/K562/WY_20210831_H3K4me3_K9Ac/4-bwCompare_bowtie2/random/H3K9ac_WT_mean.bw \
/analysis/lixf/CTBP1/ChIP/human/K562/WY_20211207_H3K27me3_K9me3/4-bwCompare_bowtie2/random/H3K9me3_WT_mean.bw \
/analysis/lixf/tracks/ATAC/GEO/rawdata/GSE213909_ATAC/4-bwCompare_bowtie2/random/ATAC_mean.bw \
/LiuLab/reference/Human/GRCh38/k100.Umap.MultiTrackMappability.bw \
 -b 3000 -a 3000 \
 --samplesLabel  K562 SY5Y  HCT116 HepG2 MCF7 A549 K27ac K4me1 K4me3 K9ac K9me3 ATAC  Mappbility \
 -o ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine.matrix.mat.gz  

plotHeatmap -m ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine.matrix.mat.gz \
            -o ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine_bilinear_for_ATAC.pdf \
            --outFileSortedRegions ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine.peak.bed \
           --heatmapHeight 10 --heatmapWidth 1.2 \
           --interpolationMethod "bilinear" \
           --colorList "#2c4688,#475E97,#5F73A5,#5E72A5,#6D7FAD,white,#EE4951,#eb1c25" "#2c4688,#475E97,#5F73A5,#5E72A5,#6D7FAD,white,#EE4951,#eb1c25" "#2c4688,#475E97,#5F73A5,#5E72A5,#6D7FAD,white,#EE4951,#eb1c25" "#2c4688,#475E97,#5F73A5,#5E72A5,#6D7FAD,white,#EE4951,#eb1c25" "#2c4688,#475E97,#5F73A5,#5E72A5,#6D7FAD,white,#EE4951,#eb1c25" "#2c4688,#475E97,#5F73A5,#5E72A5,#6D7FAD,white,#EE4951,#eb1c25" "#2c4688,white,#f7aaac,#eb1c25"  "#2c4688,white,#f7aaac,#eb1c25" "#2c4688,white,#f7aaac,#eb1c25" "#2c4688,white,#f7aaac,#eb1c25" "#2c4688,white,#f7aaac,#eb1c25" "#2c4688,#475E97,#5F73A5,#5E72A5,#6D7FAD,white,#EE4951,#eb1c25"  "#2c4688,white,#FACACB"   \
            --colorNumber 300 \
            --startLabel "5'" \
            --xAxisLabel " " \
            --endLabel "3'"  \
            --whatToShow 'heatmap and colorbar' \
            --sortUsingSamples 1 \
            --regionsLabel "PA1" "PA2" "PA3" "PA4" "PA5" "PA6" "PA7" "PA8_17" "PB" "L1M" \
            --zMax  2 2.3 2.3 2.2 2.5 2.5 10 20 20 20 50 4  1 \
            --sortUsing mean   &



plotProfile -m ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine.matrix.mat.gz \
              -out ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine_profile_grouped.pdf \
              --outFileNameData ${dir}/4-Deeptools/FLL1_profile/${prefix}/combine_profile_grouped.txt \
              --numPlotsPerRow 8 \
              --plotWidth 10 \
              --plotHeight 10 \
              --legendLocation upper-right  &
             # --colors "#bd3106" "#d9700e" "#eebe04" "#c3d6ce" "#89a6bb" "#454b87" & 



### --colorList "#2c4688,white,#f7aaac,#f16164,#eb1c25" \
} &

wait




