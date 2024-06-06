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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee1/major_point2_L1_percent/")

TPM <- fread("./TPM/RepeatMasker_Unique.TPM.tsv") 
FLL1 <- TPM %>% 
	filter(str_detect(rowname,"^L1")) %>% 
	mutate(  tmp = str_split_fixed(rowname,"_chr",2)[,2],
					 start = str_split_fixed(tmp,"_",4)[,2] %>% as.numeric ,
					 end =  str_split_fixed(tmp,"_",4)[,3] %>% as.numeric
	) %>% 
	filter(end - start >= 5999)
melt <- FLL1 %>% select(1:7) %>% melt()
meanTPM <- melt %>% mutate(cellLine = str_remove(variable,"_[0-9]")) %>% 
	group_by(cellLine,rowname) %>% 
	summarize(meanTPM = mean(value))
dcast <- dcast(meanTPM,rowname~cellLine)
FilterDcast <- dcast[apply(dcast[,-1],1,function(x) max(x) > 0.01),]
FilterDcast
mat <- FilterDcast[,-1] %>% as.matrix()
rownames(mat) <- FilterDcast$rowname
WTmat <- mat %>% t()
WTmat2 <- WTmat
col_fun = circlize::colorRamp2(c(-5,-1,0,0.3,5), colorRampPalette(c("navy", "white", "firebrick3"))(5))
dim(WTmat2)
pdf("individual_L1_plot_TPM_TPM0.01.pdf",width = 8,height = 6.7)
ComplexHeatmap::pheatmap(WTmat2,
												 #row_split = annotation_col$FactorType,
												 scale = "column",
												 col = col_fun,
												 clustering_method = "ward.D2",
												 # color = colorRampPalette(c("#8ca3bc","white","#AA4646","#9E2C2C"))(31),
												 # breaks = c(seq(0,0.9,0.1),1,seq(1.1,2,0.1),seq(2.1,3,0.1)),
												 #scale = "column",
												 name = "L1 expression \n(zscore-transformed TPM)",
												 show_colnames = F,
												 border = F,
												 #border_color = "white",
												 cellheight = 5.5,cellwidth = 0.4,
												 fontsize = 5,
												 cluster_cols = T,cluster_rows = F)
dev.off()



