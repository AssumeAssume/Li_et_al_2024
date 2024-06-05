
cd /analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/K562_WGS_INSERTION/2-Minimap2

# samtools view  -h -F 2304 -@10 WGS_K562_reporter.sam > filtered.sam
# java -jar /LiuLab/software/picard.jar SamToFastq  I=filtered.sam FASTQ=output.fastq



ref_GFP=/LiuLab/reference/Human/GRCh38/LNC5/GFP_ont.mmi
minimap2 --sam-hit-only -t 72 -ax map-ont  $ref_GFP output.fastq > GFP.sam 
samtools view -h -F 2304 -@10 GFP.sam > filtered.GFP.sam
java -jar /LiuLab/software/picard.jar SamToFastq  I=filtered.GFP.sam FASTQ=output.GFP.fastq


ref_Genome=/LiuLab/reference/Human/GRCh38/minimap2/GRCh38.p12.genome.mmi
minimap2 --sam-hit-only -t 72 -ax map-ont  $ref_Genome output.GFP.fastq > GFP.genome.sam 

samtools view -h -F 2304 -@10 GFP.genome.sam  > filtered.GFP.genome.sam 

samtools view -bhSu -@72 filtered.GFP.genome.sam |samtools sort -@ 10 - > ../3-Sam2Bam/filtered.GFP.genome.bam


cd /analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/K562_WGS_INSERTION/3-Sam2Bam

samtools index filtered.GFP.genome.bam
bamCoverage -b filtered.GFP.genome.bam \
                    --normalizeUsing CPM \
                    -p 70 \
                    --binSize 1 \
                    -o filtered.GFP.genome.bw \
                    --ignoreDuplicates 



samtools view -h filtered.GFP.genome.bam |bamToBed  -cigar > gfp.genome.bed 

bedtools intersect -a gfp.genome.bed -b /LiuLab/reference/Human/GRCh38/TE/TE_bed/LINE-1.strand.bed -f 0.9 -v  > gfp.genome.filterL1.bed

bedtools genomecov -i  gfp.genome.filterL1.bed -g /LiuLab/reference/Human/GRCh38/GRCh38.p12.genome.chrom.sizes -bg > gfp.genome.filterL1.bedgraph

bedtools cluster -i gfp.genome.filterL1.bedgraph -d 10 > gfp.genome.filterL1.cluster.bedgraph
cat gfp.genome.filterL1.cluster.bedgraph|datamash -s -g 5 max 4 |awk 'BEGIN{FS=OFS="\t"}{print $2,$1}'> max.tmp
fgrep -f max.tmp -w   gfp.genome.filterL1.cluster.bedgraph > gfp.genome.filterL1.cluster.max.bedgraph
cat gfp.genome.filterL1.cluster.max.bedgraph|awk 'BEGIN{FS=OFS="\t"}{print $0,$3-$2}'|datamash -s -g 5 min 6 -f |awk 'BEGIN{FS=OFS="\t"}NF{NF--};1' > gfp.genome.filterL1.cluster.max.L1.final.bedGarph


cat gfp.genome.filterL1.cluster.max.L1.final.bedGarph >tmp
bedtools intersect -a gfp.genome.bed -b tmp -wa -wb | sed -E -e 's/([0-9]+S).*[0-9]+[MIDNHPX=]([0-9]+S)/Left\1Right\2/g' -e 's/([0-9]+S).*[0-9]+[MIDNHPX=]/Left\1/g' -e 's/[0-9]+[MIDNHPX=].*[0-9]+[MIDNHPX=]([0-9]+S)/Right\1/g'  > tmp2

Rscript /analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/K562_WGS_INSERTION/sh_file/get_1bp_insertion.R