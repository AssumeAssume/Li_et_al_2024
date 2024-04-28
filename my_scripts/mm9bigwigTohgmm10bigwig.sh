#!/bin/bash 

bigwig=$1

bigWigToBedGraph  $bigwig ${bigwig}_tmp

liftOver ${bigwig}_tmp  /LiuLab/reference/Mouse/GRCm38/liftover_chain/mm9ToMm10.over.chain  ${bigwig}_tmp1 unMapped

cat ${bigwig}_tmp1|grep -Pi -w  "^chr[0-9XYM]+" >${bigwig}_tmp2

cat ${bigwig}_tmp2|awk 'BEGIN{FS=OFS="\t"}{getline a;$3=$3-1;print $0"\n"a}' >${bigwig}_tmp3 
bedSort ${bigwig}_tmp3 ${bigwig}_tmp4
bedtools merge -i ${bigwig}_tmp4  -c 4 -d 0 -o mean >${bigwig}_tmp5

/LiuLab/software/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${bigwig}_tmp5 /LiuLab/reference/Mouse/GRCm38/GRCm38.p6.genome.chrom.sizes ${bigwig}_mm10.bigwig

#rm -rf ./${bigwig}_tmp*rm -rf *mm10*



