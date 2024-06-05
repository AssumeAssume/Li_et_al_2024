

dir=/analysis/lixf/sgRNA/human/LINE_screening_combine/figure4_CTBP1/absolute_distance

#sizeList=(0 1 2 4 8 16 32 64 128 256 512 1024 2048 )
GeneBed=/analysis/lixf/RNA/human/K562/KO/CTBP1/5-DESeq2/KnownGene/CTBP1_upgene_removeContig.tss.sort.bed

##### K27 increased L1

queryBed=/analysis/lixf/sgRNA/human/LINE_screening_combine/figure4_CTBP1/CTBP1KO_K27ac_increase_AllL1.sort.bed
chromsize=/LiuLab/reference/Human/GRCh38/GRCh38.p12.genome.chrom.sizes

prefix=K27ac_increased_L1
outdir=${dir}/${prefix}
mkdir -p ${outdir}


bedtools closest -a ${queryBed} -b ${GeneBed} -g ${chromsize} -t first -d  >${outdir}/distance.bed





##### background expectation


chromsize=/LiuLab/reference/Human/GRCh38/GRCh38.p12.genome.chrom.sizes

prefix=background_exp
outdir=${dir}/${prefix}
mkdir -p ${outdir}

L1=/LiuLab/reference/Human/GRCh38/TE/TE_bed/LINE-1.strand.bed
bedtools intersect -a $L1 -b ${queryBed} -f 1 -r -v > ${outdir}/background.L1.tmp.bed
N=`cat ${queryBed}|wc -l`

for (( i = 1; i <= 1000; i++ )); do
	{
		shuf -n ${N}  ${outdir}/background.L1.tmp.bed |bedtools sort -g ${chromsize} >${outdir}/background.shuf.${i}.bed

	backGround=${outdir}/background.shuf.${i}.bed


	bedtools closest -a ${backGround} -b ${GeneBed} -g ${chromsize} -t first -d  >${outdir}/distance.${i}.bed
} &
done
