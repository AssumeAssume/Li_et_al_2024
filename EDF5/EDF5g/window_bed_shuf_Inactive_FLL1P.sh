dir=/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/Figure2h_RNA_K27ac
scripts=7-windowBed.sh



N=`cat /analysis/lixf/RNA/human/K562/KO/CTBP1/WY_20210306/6-LINE1M/up_p.0.05.L1.bed|wc -l`


jobs=10
[ -e /tmp/fd1 ] || mkfifo /tmp/fd1 
exec 3<>/tmp/fd1   
rm -rf /tmp/fd1 

for ((i=1;i<=${jobs};i++))
do
    echo >&3
done
for (( i = 1; i <= 100; i++ )); do
	read -u3
	{
	shuf -n ${N}  /analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/Figure2h_RNA_K27ac/Inactive_FLL1P.bed  >/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/Figure2h_RNA_K27ac/shuf_Inactive_FLL1P/shuf_L1_bed/background.shuf.${i}.bed

	backGround=/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/Figure2h_RNA_K27ac/shuf_Inactive_FLL1P/shuf_L1_bed/background.shuf.${i}.bed
	prefix=shuf_${i}
	${scripts} /analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/Figure2h_RNA_K27ac/shuf_Inactive_FLL1P/ ${prefix} $backGround -2 2 
	echo >&3
} &
done
wait
exec 3<&-
exec 3>&-