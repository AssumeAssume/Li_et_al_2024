#!/bin/bash 

bigwig=$1

bigWigToBedGraph  $bigwig tmp

liftOver tmp  /LiuLab/reference/Human/GRCh38/liftover_chain/hg19ToHg38.over.chain  tmp1 unMapped

cat tmp1|grep -Pi -w  "^chr[0-9XYM]+" >tmp2

cat tmp2|awk 'BEGIN{FS=OFS="\t"}{getline a;$3=$3-1;print $0"\n"a}' >tmp3 
bedSort tmp3 tmp4
bedtools merge -i tmp4  -c 4 -d 0 -o mean >tmp5

/LiuLab/software/kentUtils/bin/linux.x86_64/bedGraphToBigWig tmp5 /LiuLab/lixf/natureliu/KO_WT_INPUT/repeat_masker/hg38.genome.size ${bigwig/.bigwig/.hg38.bigwig}



