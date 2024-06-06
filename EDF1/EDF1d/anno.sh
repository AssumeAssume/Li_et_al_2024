

cd /analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/K562_WGS_INSERTION/4-annoPlot

homerAnno=/LiuLab/reference/Human/GRCh38/HOMER_REFERENCE/GENCODEv31.anno.final.txt
annotatePeaks.pl resolution_1bp.L1insertion.removeContig_chrY.bed hg38 -ann $homerAnno  >resolution_1bp.L1insertion.removeContig_chrY.anno.txt