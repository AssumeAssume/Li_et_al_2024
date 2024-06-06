
# cat /analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/multiple_mapped_reads_individual/L1HS_L1PA_FL_remove2L1PA13_oneL1PA8A.gtf|cut -f 9|sed 's/;.*//g'|sed 's/.*\ "//g'|sed 's/"//g' > tmp
# paste /analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/multiple_mapped_reads_individual/L1HS_L1PA_FL_remove2L1PA13_oneL1PA8A.gtf tmp|awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,"unique_id \""$10"_"$1"_"$4"_"$5"_"$7"\";"}' > L1HS_L1PA_FL_remove2L1PA13_oneL1PA8A_unique.gtf

RepeatMasker=/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/multiple_mapped_reads_individual/L1HS_L1PA_FL_unique.gtf
bam=/analysis/lixf/RNA/human/K562/KO/CTBP1/WY_20210306/1-fastq2sam2bam/bam/CTBP1_merge.sorted.bam
featureCounts -T 64 -p -t exon -g unique_id -Q 10 -f -a $RepeatMasker -o ./unique_count.featureCounts.txt ${bam} 
awk '{if($0 !~ "^#"){print $1"\t"$7}}' ./unique_count.featureCounts.txt > ./unique_count.rawcounts.tsv

featureCounts -T 64 -p -t exon -g unique_id -Q 0 -f -M --fraction -a $RepeatMasker -o ./multiple_count.featureCounts.txt ${bam}
awk '{if($0 !~ "^#"){print $1"\t"$7}}' ./multiple_count.featureCounts.txt > ./multiple_count.rawcounts.tsv




