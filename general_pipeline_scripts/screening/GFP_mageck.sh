#!/bin/bash
conda activate mageck-vispr
# mageck-vispr init /analysis/lixf/transmit/esc-testdata/GFP --reads /LiuLab/rawdata/N2009231_YXH_80-542388554_SEQ/201103-X4A_L006/1-plus_combined_R1.fastq.gz /LiuLab/rawdata/N2009231_YXH_80-542388554_SEQ/201103-X4A_L006/2-plus_combined_R1.fastq.gz /LiuLab/rawdata/N2009231_YXH_80-542388554_SEQ/201103-X4A_L006/1-NC_combined_R1.fastq.gz /LiuLab/rawdata/N2009231_YXH_80-542388554_SEQ/201103-X4A_L006/2-NC_combined_R1.fastq.gz


cd /analysis/lixf/transmit/esc-testdata/GFP

# snakemake -n

# snakemake --cores 30 &

cat   ./results/count/all.count_normalized.txt|sed '1d'|datamash -s -g 2 sum 3,4,5,6  |sort -k 1|sed '1igene\t1-treat\t2-treat\t1-control\t2-control' > ./results/count/sum_normalized_count.tsv


# cp ./results/test/mle.gene_summary.txt  ./results/test/MaGeCK_L1_GFP_positive.gene_summary.txt

#vispr server results/*.vispr.yaml --host 10.10.31.12

cd /analysis/lixf/transmit/esc-testdata/GFP_2

 snakemake -n

 snakemake --cores 30 &
 cp ./results/test/mle.gene_summary.txt  ./results/test/MaGeCK_L1_GFP_negative.gene_summary.txt

#vispr server results/*.vispr.yaml --host 10.10.31.12
cat   ./results/count/all.count_normalized.txt|sed '1d'|datamash -s -g 2 sum 3,4,5,6  |sort -k 1|sed '1igene\t1-treat\t2-treat\t1-control\t2-control' > ./results/count/sum_normalized_count.tsv


