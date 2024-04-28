#!/bin/bash

# THREAD_NUM for a time 
THREAD_NUM=1

#定义描述符为9的FIFO管道
[ -e /tmp/fifo9 ] || mkfifo /tmp/fifo9 
exec 9<>/tmp/fifo9   
rm -rf /tmp/fifo9 

#预先写入指定数量的空格符，一个空格符代表一个进程
for ((i=0;i<$THREAD_NUM;i++))
do
    echo >&9
done

export PATH="/analysis2/software/miniconda3/envs/tf/bin:$PATH" 
export LD_LIBRARY_PATH="/analysis2/software/miniconda3/envs/tf/lib:$LD_LIBRARY_PATH"
export TMPDIR=/analysis2/lixf/Hi-C/human//scripts/tmpdir

cd /analysis2/lixf/Hi-C/human/

mkdir -p 2-DeepLoop

DeepLoop=/analysis2/software/DeepLoop-master
HiCorr=/analysis2/software/HiCorr-master

hg19=/analysis2/reference/Human/hg19/bowtie2/hg19
hg19Ind=/analysis2/reference/Human/hg19/hg19.fa.fai
HiCorrPath=/analysis2/software/HiCorr-master
lib=/analysis2/software/HiCorr-master/lib
bed=/analysis2/software/HiCorr-master/ref/DPNII/hg19.DPNII.frag.bed # <chr> <start> <end> <frag_id>
anchorDir=${DeepLoop}/DeepLoop_models/ref/hg19_DPNII_anchor_bed/

threads=50
dir=/analysis2/lixf/Hi-C/human///2-DeepLoop

metadata=/analysis2/lixf/Hi-C/human///hic_metadata_alias.txt
declare -A array
data_grep=$
eval $(cat ${metadata} |grep -Pi "${data_grep}"|sed '/^$/d'|datamash -s -g 1 collapse 2 |awk -F "\t" '{print "array["$1"]="$2}')

for expt in ${!array[@]}; do
  read -u9
    {
        ##### Pre HiCorr data processing
        mkdir -p ${dir}/${expt}/0-PreHiCorr
        mkdir -p ${dir}/${expt}/1-HiCorr
        mkdir -p ${dir}/${expt}/2-DeepLoop
        mkdir -p ${dir}/${expt}/3-coolerFormat/heatmaps

        outdir=${dir}/${expt}/0-PreHiCorr
      ##  JuicerfastqDir=/analysis2/lixf/Hi-C/human///1-Juicer
        readsarray=(`echo ${array[$expt]}|tr ',' ' '`)
        Count=$(echo "${#readsarray[@]}" | bc)

        for ((k=1;k<=$Count;k++));do

            fq1=${array[$expt]}_1.fq.gz
            fq2=${array[$expt]}_2.fq.gz

            # #### run mapping with bowtie (hg19 build)
            bowtie2 --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end  --rg-id BMG --rg SM:${expt}_run${k}_R1 -p ${threads} -x $hg19 -U $fq1 |  samtools view -bhSu  | samtools sort -@ ${threads} -m 4G -O bam -T $expt.run${k}.R1 -o ${outdir}/$expt.run${k}.R1.sorted.bam &
            bowtie2 --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end  --rg-id BMG --rg SM:${expt}_run${k}_R2 -p ${threads} -x $hg19 -U $fq2 |  samtools view -bhSu  | samtools sort -@ ${threads} -m 4G -O bam -T $expt.run${k}.R2 -o ${outdir}/$expt.run${k}.R2.sorted.bam &

            samtools sort -@ ${threads} -n -T $expt.run${k}.R1 -o ${outdir}/$expt.run${k}.R1.name.sorted.bam ${outdir}/$expt.run${k}.R1.sorted.bam &
            samtools sort -@ ${threads} -n -T $expt.run${k}.R2 -o ${outdir}/$expt.run${k}.R2.name.sorted.bam ${outdir}/$expt.run${k}.R2.sorted.bam &
            wait

            samtools index  ${outdir}/$expt.run${k}.R1.name.sorted.bam  &
            samtools index  ${outdir}/$expt.run${k}.R2.name.sorted.bam &
            
             wait


    

             python $lib/mergeSAM.py -q 0 -t -v -f ${outdir}/$expt.run${k}.R1.name.sorted.bam -r ${outdir}/$expt.run${k}.R2.name.sorted.bam -o ${outdir}/$expt.run${k}.merge.bam

        done &
        wait

        ### merge files belong to the same biological replicate and remove duplicates
        samtools merge ${outdir}/$expt.bam ${outdir}/$expt.run*.merge.bam


        samtools sort -@ ${threads} -T ${expt} ${outdir}/$expt.bam | samtools view - | perl $lib/remove_dup_PE_SAM_sorted.pl | samtools view -bS -t $hg19Ind -o - - > ${outdir}/$expt.sorted.nodup.bam

        ### map reads pair to fragment pairs
        cd ${dir}/${expt}/0-PreHiCorr
        ${HiCorr}/HiCorr Bam-process-DpNII ${outdir}/$expt.sorted.nodup.bam ${expt} 150 hg19 DPNII
        
        ### Run HiCorr
        cd ${dir}/${expt}/1-HiCorr
        ${HiCorr}/HiCorr DPNII ${dir}/${expt}/0-PreHiCorr/${expt}.cis.frag_loop ${dir}/${expt}/0-PreHiCorr/${expt}.trans.frag_loop ${expt} hg19


        ## Run DeepLoop
        for chr in {1..22..1} X; do
            python3 ${DeepLoop}/prediction/predict_chromosome.py --full_matrix_dir ${dir}/${expt}/1-HiCorr/HiCorr_output \
                                                      --input_name anchor_2_anchor.loop.chr${chr} \
                                                      --h5_file ${DeepLoop}/DeepLoop_models/CPGZ_trained/2.4M.h5 \
                                                      --out_dir ${dir}/${expt}/2-DeepLoop/ \
                                                      --anchor_dir ${anchorDir} \
                                                      --chromosome chr${chr} \
                                                      --small_matrix_size 128 \
                                                      --step_size 128 \
                                                      --dummy 5 \
                                                      --val_cols obs exp
        done


        ## convert to cooler
        for chr in {1..22..1} X; do
            python3 ${DeepLoop}/utils/convert_to_cooler.py --anchor_dir ${anchorDir} \
                                                      --loop_dir ${dir}/${expt}/2-DeepLoop/ \
                                                      --out_file ${dir}/${expt}/3-coolerFormat/${expt}.chr${chr}.cool \
                                                      --col_names a1 a2 enhance \
                                                      --cooler_col enhance \
                                                      --single_chrom chr${chr}
        done

        python3 ${DeepLoop}/utils/convert_to_cooler.py --anchor_dir ${anchorDir} \
                                                      --loop_dir ${dir}/${expt}/2-DeepLoop/ \
                                                      --out_file ${dir}/${expt}/3-coolerFormat/${expt}.enhance.cool \
                                                      --col_names a1 a2 enhance \
                                                      --cooler_col enhance &

        cd ${dir}/${expt}/3-coolerFormat
        python3 ${DeepLoop}/utils/convert_to_cooler.py --anchor_dir ${anchorDir} \
                                                      --loop_dir ${dir}/${expt}/2-DeepLoop/ \
                                                      --out_file ${dir}/${expt}/3-coolerFormat/${expt}.enhance.5kb.cool \
                                                      --col_names a1 a2 enhance \
                                                      --cooler_col enhance \
                                                      --bin_size 5000 \
                                                      --min_val 1.0 \
                                                      --force_bin_size &

        python3 ${DeepLoop}/utils/convert_to_cooler.py --anchor_dir ${anchorDir} \
                                                      --loop_dir ${dir}/${expt}/2-DeepLoop/ \
                                                      --out_file ${dir}/${expt}/3-coolerFormat/${expt}.enhance.10kb.cool \
                                                      --col_names a1 a2 enhance \
                                                      --cooler_col enhance \
                                                      --bin_size 10000 \
                                                      --min_val 1.0 \
                                                      --force_bin_size &


    echo >&9
    } &
done
wait

exec 9<&-
exec 9>&-

