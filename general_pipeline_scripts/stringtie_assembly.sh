

ref="/LiuLab/reference/Human/GRCh38/GTF/gencode.v31.chr_patch_hapl_scaff.annotation.gff"
dir=/analysis/lixf/CRISPRa_i/NCCIT/RNA/fusion_transcript

mkdir -p ${dir}/1-stringtie_output
mkdir -p ${dir}/2-taco_output
# for bam in `ls ${dir}/0-bamlist|grep .sorted.bam$`
# do
# 	stringtie ${dir}/0-bamlist/${bam} -p 8 -j 2 -s 5 -f 0.05 -c 2 -G ${ref} -o ${dir}/1-stringtie_output/`basename ${bam} .sorted.bam`_stringtie_output.gtf &
# done
# wait

# readlink -f ${dir}/1-stringtie_output/*gtf > ${dir}/2-taco_output/gtf_to_merge.txt

# stringtie --merge -p 30 -G ${ref} -i -o ${dir}/2-taco_output//stringtie_merged.gtf  ${dir}/2-taco_output/gtf_to_merge.txt

taco_run -o ${dir}/2-taco_output/merge_gtf -p 40 ${dir}/2-taco_output/gtf_to_merge.txt

taco_refcomp -o ${dir}/2-taco_output/anno -p 40 -r /LiuLab/reference/Human/GRCh38/GTF/gencode.v31.chr_patch_hapl_scaff.annotation.gtf -t ${dir}/2-taco_output/merge_gtf/assembly.gtf