#!/bin/bash


RepeatMasker=/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/second_rebuttle/Unique_Multiple_percent/L1HS_L1PA_FL.gtf
bam=/analysis/lixf/RNA/human/K562/KO/CTBP1/WY_20210306/1-fastq2sam2bam/bam/CTBP1_merge.sorted.bam
featureCounts -T 64 -p -t exon -g gene_id -Q 10  -a $RepeatMasker -o ./unique_count.featureCounts.txt ${bam} 
awk '{if($0 !~ "^#"){print $1"\t"$7}}' ./unique_count.featureCounts.txt > ./unique_count.rawcounts.tsv

featureCounts -T 64 -p -t exon -g gene_id -Q 0 -M --fraction -a $RepeatMasker -o ./multiple_count.featureCounts.txt ${bam}
awk '{if($0 !~ "^#"){print $1"\t"$7}}' ./multiple_count.featureCounts.txt > ./multiple_count.rawcounts.tsv




