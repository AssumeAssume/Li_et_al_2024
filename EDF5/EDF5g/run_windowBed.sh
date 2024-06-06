
dir=/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/Figure2h_RNA_K27ac
scripts=7-windowBed.sh

bedList=(/analysis/lixf/RNA/human/K562/KO/CTBP1/WY_20210306/6-LINE1M/up_p.0.05.L1.bed
	
	)

for bed in ${bedList[@]}; do
	prefix=`echo $bed|rev|cut -d '/' -f1|rev|sed 's/.bed//g'` 
	${scripts} ${dir} ${prefix} $bed -2 2 &
done


