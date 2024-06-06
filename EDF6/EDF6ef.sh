#!/bin/bash

dir=/analysis/lixf/KMT2D/ChIP/human/K562//
mkdir -p ${dir}/4-Deeptools/TE_heatmap_combine_K27ac_desc/
#########  N-dox  peak  #################
arranged_peak=/analysis/lixf/KMT2D/ChIP/human/K562//5-DESeq2/WT_KO_merge_peak/combine_peak_sort_log2FC_desc.bed
computeMatrix reference-point \
 --referencePoint center \
 -p 72 \
 --missingDataAsZero \
 --skipZeros \
 -R ${arranged_peak} \
 --binSize 500 \
 -S /analysis/lixf/KMT2D/ChIP/human/K562/4-bwCompare_bowtie2/H3K4me1_WT_mean.bw \
 /analysis/lixf/KMT2D/ChIP/human/K562/4-bwCompare_bowtie2/H3K4me1_KO_mean.bw \
 /analysis/lixf/KMT2D/ChIP/human/K562/4-bwCompare_bowtie2/H3K4me1_KO_vs_WT_subtract.bw \
/analysis/lixf/KMT2D/ChIP/human/K562/4-bwCompare_bowtie2/H3K27ac_WT_mean.bw \
 /analysis/lixf/KMT2D/ChIP/human/K562/4-bwCompare_bowtie2/H3K27ac_KO_mean.bw \
/analysis/lixf/KMT2D/ChIP/human/K562/4-bwCompare_bowtie2/H3K27ac_KO_vs_WT_subtract.bw \
      /LiuLab/reference/Human/GRCh38/TE/TE_bw/L1.FullLength.bw \
      /LiuLab/reference/Human/GRCh38/TE/TE_bw/L1M.bed.bw \
      /LiuLab/reference/Human/GRCh38/TE/TE_bw/LTR.bw \
      /LiuLab/reference/Human/GRCh38/TE/TE_bw/SINE.bw \
      /LiuLab/reference/Human/GRCh38/TE/TE_bw/DNA.bw \
 -b 15000 -a 15000 \
 --samplesLabel  H3K4me1_WT H3K4me1_KO  subtract H3K27ac_WT H3K27ac_KO K27_subtract  FLL1 L1M LTR SINE DNA  \
 -o ${dir}/4-Deeptools/TE_heatmap_combine_K27ac_desc/selected.matrix.mat.gz \
 --outFileSortedRegions ${dir}/4-Deeptools/TE_heatmap_combine_K27ac_desc/KMT2D_TE.peak.bed



plotHeatmap -m ${dir}/4-Deeptools/TE_heatmap_combine_K27ac_desc/selected.matrix.mat.gz \
            -o ${dir}/4-Deeptools/TE_heatmap_combine_K27ac_desc/KMT2D_TE_selected_FLL1_0.5_K27ac_KOwt_colorbar.pdf \
           --heatmapHeight 20 --heatmapWidth 3 \
            --colorList "#2c4688,#475E97,#5F73A5,#5E72A5,#6D7FAD,white,#EE4951,#eb1c25" "#2c4688,#475E97,#5F73A5,#5E72A5,#6D7FAD,white,#EE4951,#eb1c25"    "#2c4688,white,#eb1c25" "#2c4688,#475E97,#5F73A5,#5E72A5,#6D7FAD,white,#EE4951,#eb1c25" "#2c4688,#475E97,#5F73A5,#5E72A5,#6D7FAD,white,#EE4951,#eb1c25"    "#2c4688,white,#EF4B52,#eb1c25"   "#2c4688,#EE3E46,#ED3C43,#ED3941,#ED373F,#ED353D,#ED323B,#EC3038,#EC2E36,#EC2C34,#EC2932,#EC2730,#EB252D,#EB222B,#EB2029,#EB1E27,#EB1C25"  "#2c4688,white, #F26F75,#F16369,#F0575E,#EF4B52,#EE3F47,#ED333B,#EC2730,#EB1C25" "#2c4688,white, #F26F75,#F16369,#F0575E,#EF4B52,#EE3F47,#ED333B,#EC2730,#EB1C25" "#2c4688,white, #F26F75,#F16369,#F0575E,#EF4B52,#EE3F47,#ED333B,#EC2730,#EB1C25" "#2c4688,white, #F26F75,#F16369,#F0575E,#EF4B52,#EE3F47,#ED333B,#EC2730,#EB1C25" "#2c4688,white, #F26F75,#F16369,#F0575E,#EF4B52,#EE3F47,#ED333B,#EC2730,#EB1C25" "#2c4688,white, #F26F75,#F16369,#F0575E,#EF4B52,#EE3F47,#ED333B,#EC2730,#EB1C25" \
            --colorNumber 300 \
            --refPointLabel " " \
            --xAxisLabel " " \
            --whatToShow 'heatmap and colorbar' \
            --zMin 0 0 -10  0 0  -5 0 0 0 0 0 0 0 0 0 0  --zMax 130 130  10 50 50  10 0.5 0.5 0.5 1 0.5  \
           --sortRegions "keep" &