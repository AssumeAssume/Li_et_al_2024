library(tidyverse)
library(magrittr)
library(karyoploteR)
theme_set(theme_bw())
rm(list=ls())
filter <- dplyr::filter
select <- dplyr::select
reduce <- purrr::reduce
options(stringsAsFactors = FALSE)
set.seed(629)


setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/L1_regulated_gene_cluster")

DistanceDf <- fread("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/mouseCRISPRai/Gene_cloest_L1/Gene_6K_L1.distance.txt")
geneList2 <- fread("/analysis/lixf/CRISPRa_i/mESC/RNA/public_gene_sets/2c_genesets.txt")
ASO2c <- fread("/analysis/lixf/CRISPRa_i/mESC/RNA/ASO2c/5-DESeq2/KnownGene/Combine_removeL1_1Deseq2.L1_vs_SCR.result.txt")
fusionDistance <- DistanceDf %>% filter(V13==0)

geneBed <- fread("/LiuLab/reference/Mouse/GRCm38/GTF/gencode.Allgene.strand.bed")
tmp <- ASO2c %>% 
	filter(rowname %in% geneList2$GSE33923 ) %>% 
	mutate(group = case_when(rowname %in% fusionDistance$V4 ~ "fusion",
													 padj < 0.1 & log2FoldChange < 0 ~ "L1-regulated",
													 TRUE ~ 'others'
													 )) 
colnames(geneBed)[4] <- "rowname"
bed <- geneBed %>%
	inner_join(tmp) %>% 
	filter(str_detect(V1,"chr"))

l1_regulated_bed <- bed %>% filter(group == "L1-regulated") %>% select(1:4)
colnames(l1_regulated_bed) <- c("chr","start","end","names")
table(c(l1_regulated_bed %>% pull(names),df %>% pull(names))) %>% sort
df <- bed %>% filter(group == "others") %>% select(1:4)
colnames(df) <- c("chr","start","end","names")

l1genes <- toGRanges(l1_regulated_bed)
other2C <- toGRanges(df)
combinedf <- bed %>%  select(1:4)
colnames(combinedf) <- c("chr","start","end","names")
combineLbael <- toGRanges(combinedf)
pdf("localization_2.pdf",height=17.5,width = 15)
kp <- plotKaryotype(genome = "mm10",plot.type=1, main="L1-regulated 2C genes")
kpPlotMarkers(kp,data=combineLbael,labels=combineLbael$names, text.orientation = "vertical",
							line.color="#575252", label.color="black",label.dist = 0.002,
							r1=0.1, cex=0.5, adjust.label.position = T,max.iter=1000)

dev.off()
