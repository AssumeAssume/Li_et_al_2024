outdir=/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/intron_splicing/intron_count
bed12=/LiuLab/reference/Human/GRCh38/GTF/gencode.v31.hg38.bed12
#bamList=/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/intron_splicing/intron_count/selected.bam.txt
bamList=/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/intron_splicing/intron_count/ENCODE_redoBAM.txt

mkdir -p ${outdir}/read_distribution
declare -A array

eval $(cat ${bamList}|sed '/^$/d'|awk -F "\t" '{print "array["$1"]="$2}')

thread=52
jobs=10
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

read_distribution.py  -i ${array[$key]} -r ${bed12}  >  ${outdir}/read_distribution/${key}.txt

echo >&3
             } & 



done
wait
exec 3<&-
exec 3>&-