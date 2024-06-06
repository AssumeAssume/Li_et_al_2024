outdir=/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/intron_splicing/SPLICEq

#bamList=/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/intron_splicing/intron_count/selected.bam.txt
bamList=/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/intron_splicing/SPLICEq/redo.bamList.txt
anno=/LiuLab/reference/Human/GRCh38/GTF/gencode.v31.chr_patch_hapl_scaff.annotation.gtf


mkdir -p ${outdir}/SPLICE_efficiency/
mkdir -p ${outdir}/SE_stat/
declare -A array

eval $(cat ${bamList}|sed '/^$/d'|awk -F "\t" '{print "array["$1"]="$2}')

thread=20
jobs=14
[ -e /tmp/fd1 ] || mkfifo /tmp/fd1 
exec 3<>/tmp/fd1   
rm -rf /tmp/fd1 

for ((i=1;i<=${jobs};i++))
do
    echo >&3
done
# clean data

for key in ${!array[@]}
    do
read -u3
        {


####### bam2bigwig
/LiuLab/software/miniconda3/bin/SPLICE-q.py -b ${array[$key]} -g ${anno} -o ${outdir}/SPLICE_efficiency/${key}.tsv -p ${thread}

cat ${outdir}/SPLICE_efficiency/${key}.tsv|datamash --headers mean 14 median 14 q1 14 q3 14  > ${outdir}/SE_stat/${key}.stat.tsv

# featureCounts -T ${thread} -p -t exon -g unique_id -f  -Q 0 -M \
#                -a $RepeatMaskerUnique  -o ${outdir}/3-featureCounts/RepeatMasker_Unique_Multimap/${key}.L1_combine.txt ${array[$key]}
# awk '{if($0 !~ "^#"){print $1"\t"$7}}' ${outdir}/3-featureCounts/RepeatMasker_Unique_Multimap/${key}.L1_combine.txt > ${outdir}/4-RawCount/RepeatMasker_Unique_Multimap/${key}.L1_combine.tsv


# featureCounts -T ${thread} -p -t exon -g unique_id -f  -Q 0 -M --fraction \
#                -a $RepeatMaskerUnique  -o ${outdir}/3-featureCounts/RepeatMasker_Unique_Multimap_Fraction/${key}.L1_combine.txt ${array[$key]}
# awk '{if($0 !~ "^#"){print $1"\t"$7}}' ${outdir}/3-featureCounts/RepeatMasker_Unique_Multimap_Fraction/${key}.L1_combine.txt > ${outdir}/4-RawCount/RepeatMasker_Unique_Multimap_Fraction/${key}.L1_combine.tsv

echo >&3
             } & 



done
wait
exec 3<&-
exec 3>&-





