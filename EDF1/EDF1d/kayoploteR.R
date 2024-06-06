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


setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/K562_WGS_INSERTION/4-annoPlot")


position <- fread("./resolution_1bp.L1insertion.removeContig_chrY.bed",header = F) %>% 
	filter(V4>3)


df <- as.data.frame(position %>% select(1:3))
names(df) <- c("chr","start","end")
pdf("localization_Coverage_3.pdf",height=3.5,width = 4)
kp <- plotKaryotype(genome = "hg38",plot.type=1, main="LINE1 reporter insertion site")

kpPlotRegions(kp, data=df,border = "#FFB318")

dev.off()
