

dir=/analysis/lixf/tracks/STARR/several_cellLine

###### ls ${dir}/3-MACS2_bowtie2|tr '\n' ' '
outdir=${dir}/11-annoEnrichment_family_subfamily


mkdir -p ${outdir}/subfamily_enrichment
mkdir -p ${outdir}/family_enrichment
homerAnno=/LiuLab/reference/Human/GRCh38/HOMER_REFERENCE/GENCODEv31.anno.final.txt
for bed in `ls ${dir}/merged_bed`
do

		{
		###  peaks overlap with  LINE-1
		prefix=${bed%%_*}

		annotatePeaks.pl ${dir}/merged_bed/${bed} \
		hg38 \
		-ann $homerAnno  \
		-annStats ${outdir}/subfamily_enrichment/${prefix}_enrichment.ann.txt > ${outdir}/subfamily_enrichment/${prefix}_annotation.txt

		annotatePeaks_family.pl ${dir}/merged_bed/${bed} \
		hg38 \
		-ann $homerAnno  \
		-annStats ${outdir}/family_enrichment/${prefix}_enrichment.ann.txt > ${outdir}/family_enrichment/${prefix}_annotation.txt

	} &
	# bedtools igv -path /analysis/lixf/KDM2B/ChIP/human/NCCIT/YXH_20210420/igv/batch/SF_C_L1P -sess /analysis/lixf/IGV/KDM2B/NCCIT_track_without_histon.xml -slop 10000 -i ./intersect_L1P.bed >../../igv/batch/SF_C_L1P/batch.sh

	# bedtools igv -path /analysis/lixf/KDM2B/ChIP/human/NCCIT/YXH_20210420/igv/batch/SF_C_L1P_DE -sess /analysis/lixf/IGV/KDM2B/NCCIT_track_without_histon.xml -slop 10000 -i ./DE_binding_L1P.bed >../../igv/batch/SF_C_L1P_DE/batch.sh
	
done
wait
