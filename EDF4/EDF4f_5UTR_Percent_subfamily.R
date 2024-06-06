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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee1/major_point4_STARR")
LINE1factor <- fread("/LiuLab/reference/Human/GRCh38/TE/UCSCrmsk/L1.div.txt",header = F) 

starr_FL <- fread("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee1/major_point4_STARR/norm_FL_include_percent.tsv") %>% dft()
fl_name <- rownames(starr_FL)


starr_nonFL <- fread("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee1/major_point4_STARR/norm_nonFL_include_percent.tsv") %>% dft()
nonfl_name <- rownames(starr_nonFL)

L1 <- fread("./5UTR.L1.true.bed") %>% rename(subfamily = V4)
L1_cn <- fread("/LiuLab/reference/Human/GRCh38/TE/TE_bed/L1_copynumber_df.tsv") %>% select(-4)
colnames(L1_cn) <- c("subfamily","total","FL")
meltCN <- L1_cn %>% mutate(`non-FL`=total-FL) %>% select(subfamily,FL,`non-FL`) %>%  melt()
colnames(meltCN)[2] <- "group"
test <- L1 %>% 
	mutate(group = ifelse(V3-V2 >=6000,"FL","non-FL")) %>% 
	group_by(group,subfamily) %>% 
	count() %>% ungroup %>% full_join(meltCN)
plotdf <- test %>% 
	mutate(percent = n/value,
				 subfamily = factor(subfamily,levels = LINE1factor$V1))


tmp <- plotdf %>% 
	filter(group== "FL" & subfamily %in% fl_name) %>% 
	dcast(subfamily~group,value.var = "percent")
tmp[is.na(tmp)] <- 0
mat <- tmp[,-1] %>% as.data.frame() 
rownames(mat) <- tmp$subfamily
col_fun = circlize::colorRamp2(c(0,0.03,0.07,0.09,0.1,0.13,1), colorRampPalette(rev(brewer.pal(n = 7, name =	"RdYlGn")))(7))

annotation_col <- data.frame(Group = rep(c("FL"),each=nrow(mat)))
label_L1 <- c("L1HS","L1PA2","L1PA3","L1PA4")
L1_position <- match(label_L1,rownames(mat))
rowanno = rowAnnotation(foo = anno_mark(at = L1_position, 
																				labels = label_L1,labels_gp = gpar(fontsize = 6)))
p1 <- ComplexHeatmap::pheatmap(mat,
															 col = col_fun,
															 border =F,cellheight = 1,cellwidth = 5.5,
															 right_annotation = rowanno,
															 annotation_row = annotation_col,
															 cluster_cols = F,fontsize = 5,show_rownames = F,
															 cluster_rows = F)

tmp <- plotdf %>% 
	filter(group== "non-FL" & subfamily %in% nonfl_name) %>% 
	dcast(subfamily~group,value.var = "percent")
tmp[is.na(tmp)] <- 0
mat <- tmp[,-1] %>% as.data.frame()
rownames(mat) <- tmp$subfamily
col_fun = circlize::colorRamp2(c(0,0.15,0.2,0.3,1), colorRampPalette(rev(brewer.pal(n = 7, name =	"RdYlGn")))(5))
annotation_col <- data.frame(Group = rep(c("Non-FL"),each=nrow(mat)))
label_L1 <- c("L1PBa1","L1P2","L1PBa","L1P3b")
L1_position <- match(label_L1,rownames(mat))
rowanno = rowAnnotation(foo = anno_mark(at = L1_position, 
																				labels = label_L1,labels_gp = gpar(fontsize = 6)))
p2 <- ComplexHeatmap::pheatmap((mat),
															 annotation_row = annotation_col,
															 col = col_fun,cellheight = 1,cellwidth = 5.5,
															 border = F,fontsize = 5,
															 right_annotation = rowanno,
															 cluster_cols = F,show_rownames = F,
															 cluster_rows = F)
p2
p1%v%p2
pdf("5UTR_percent_selected.pdf",width = 3,height = 5)
p <- p1%v%p2
draw(p,ht_gap = unit(1, "mm"))
dev.off()
