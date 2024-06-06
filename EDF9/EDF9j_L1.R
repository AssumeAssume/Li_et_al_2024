library(tidyverse)
library(magrittr)
library(data.table)
library(glue)
theme_set(theme_bw())
rm(list=ls())
filter <- dplyr::filter
select <- dplyr::select
desc <- dplyr::desc
reduce <- purrr::reduce
rename <- dplyr::rename
count <- dplyr::count
readr::local_edition(1)
options(stringsAsFactors = FALSE)
set.seed(629)
source("/analysis/lixf/my_script/my_R_function.R")





setwd("/analysis2/lixf/Hi-C/mouse/mESC/CRISPRai/BLY_20230926//11-coolpup_juicerLoop_CombineRep/FLL1_5K_CRISPRi_NCloop/")


CRISPRi <- fread("./Unbanlanced_Loop.CRISPRi.loop.txt",skip = "35") %>% as.matrix()
iNC <- fread("./Unbanlanced_Loop.iNC.loop.txt",skip = "35")%>% as.matrix()

library(ComplexHeatmap)
colorRampPalette(rev(brewer.pal(n = 7, name =   "RdBu")))(3)
col_fun = circlize::colorRamp2(c(0.5,1,2.0,3.7,4.5), colorRampPalette(rev(brewer.pal(n = 7, name =	"RdBu")))(5))
hCRISPRi <- ComplexHeatmap::pheatmap(CRISPRi,cluster_rows = F,cluster_cols = F,
															 col = col_fun,
															 border=F,cellwidth = 5,cellheight = 5,
															 #main="2C",
															 show_rownames = F,show_colnames = F)

hNC <- ComplexHeatmap::pheatmap(iNC,cluster_rows = F,cluster_cols = F,
																col = col_fun,
																border=F,cellwidth = 5,cellheight = 5,
																# main="8C",
																show_rownames = F,show_colnames = F)


pdf("NC_CRISPRi_heatmap.pdf",width = 5,height = 3)
hNC+hCRISPRi
dev.off()
	