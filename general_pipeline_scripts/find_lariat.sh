
##### There must be gemnome.fa file and ss_table file under the index dir

######### Generate ss_table from ncbi_refseq 
#### cat /analysis2/reference/Human/GRCh38/GTF/hg38_ncbiRefseq.txt|awk 'BEGIN{FS=OFS="\t"}{id=$2; chr=$3; txstart=$5; txend=$6;exonstart=$10; exonend=$11; strand=$4; exonCount=$9; print id,chr,strand,txstart,txend,exonCount,exonstart,exonend}' > /analysis2/reference/Human/GRCh38/bowtie/GRCh38.p12.genome_ss_table.txt

dir=/analysis/lixf/RNA/human/K562/KO/DBR1/WY_20210306/
### single_end ; paired
layout=paired

index=/LiuLab/reference/Human/GRCh38/bowtie/GRCh38.p12.genome

mkdir -p ${dir}/6-find_larit

outdir=${dir}/6-find_larit/


for i in `ls ${dir}/0-clean_data/${layout}/|grep gz$`; do
	#statements
{
	fastq=${dir}/0-clean_data/${layout}/${i}
	 echo -e "########################## \n analysising $fastq \n##############################"
	

	if [[ ! -e ${fastq/.gz/}  ]]; then
		#statements
		pigz -d -k -p 8 $fastq
	fi

	prefix=`echo ${fastq}|rev| cut -d "/" -f  1|rev| sed 's/.fq.gz//'`

	perl /analysis2/lixf/my_script/findlariats/lariat_scripts/find_lariats.pl -f ${fastq/.gz/} -i $index -o ${outdir}/${prefix} -l 151 -m 8

} &

done

wait

### normalize larait reads

cd  ${dir}/0-QC/flag_stat


for i in `ls|grep stat$`
do
awk 'NR==1{total=$1}NR==5{mapped=$1;gsub("\\(", "",$5);rate=$5;gsub(".sorted.stat","",FILENAME); print total,mapped,rate,FILENAME}' $i >>summary.txt
done

cd ${dir}/6-find_larit

cat */lariat_data_table.txt > combine_larit_table.txt



