cd /analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/STARR


cCRE_gtf=/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/STARR/cCRE_randomsequence/combined_cCRE_randomSequence.gtf


mkdir -p ./RS_cCRE_subfamily
mkdir -p ./RS_cCRE_individual

for i in `cat bamList`; do
	prefix=`basename $i .sorted.rmdup.bam`


	featureCounts -T 64 -p -t exon -g unique_id -Q 10 -f -a $cCRE_gtf -o ./RS_cCRE_individual/${i}.5UTR.txt /analysis/lixf/tracks/STARR/several_cellLine/rawdata_analysis/1-bowtie2/bam/random/${i}
	awk '{if($0 !~ "^#"){print $1"\t"$7}}' ./RS_cCRE_individual/${i}.5UTR.txt > ./RS_cCRE_individual/${i}.rawcounts.5UTR.tsv

	featureCounts -T 64 -p -t exon -g gene_id -Q 10  -a $cCRE_gtf -o ./RS_cCRE_subfamily/${i}.5UTR.txt /analysis/lixf/tracks/STARR/several_cellLine/rawdata_analysis/1-bowtie2/bam/random/${i}
	awk '{if($0 !~ "^#"){print $1"\t"$7}}' ./RS_cCRE_subfamily/${i}.5UTR.txt > ./RS_cCRE_subfamily/${i}.rawcounts.5UTR.tsv

done




