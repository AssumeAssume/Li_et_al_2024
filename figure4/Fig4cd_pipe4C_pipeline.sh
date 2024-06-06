dir=/analysis/lixf/CRISPRa_i/NCCIT/4C/

pipe4C=/analysis/lixf/my_script/pipe4C/pipe4C.R
configFile=/analysis/lixf/my_script/pipe4C/conf.yml

mkdir -p ${dir}/pipe4C_bowtie2
Rscript ${pipe4C} --vpFile=${dir}/scripts/VPinfo.txt --fqFolder=${dir}/0-rawdata --outFolder=${dir}/pipe4C_bowtie2_keepAll --confFile=${configFile} --cores 50 --plot --wig --genomePlot  --tsv  --bins
