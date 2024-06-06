library(tidyverse)
library(magrittr)
library(data.table)
library(glue)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
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


LINE1factor2 <- fread("/LiuLab/reference/Human/GRCh38/TE/LINE_1_factor.txt",header = F) 

LINE1factor <- fread("/LiuLab/reference/Human/GRCh38/TE/UCSCrmsk/L1.div.txt",header = F) 

file <- dir("/analysis/lixf/tracks/STARR/several_cellLine/overlap_L1","overlap_L1.bed",full.names = T)
names(file) <- str_split_fixed(str_remove_all(file,".*\\/"),"_",2)[,1]
file
df_combine <- imap(file,~fread(.x) %>% mutate(cellLine = .y)) %>% reduce(rbind) %>%  mutate(length = V3-V2) 

L1_cn <- fread("/LiuLab/reference/Human/GRCh38/TE/TE_bed/L1_copynumber_df.tsv") %>% select(-4)
colnames(L1_cn) <- c("subfamily","total","FL")
#fwrite(df_combine,"L1_starr_overlap.bed",sep = '\t')
L1_cn_new <- L1_cn %>% mutate(`Non_FL` = total -  FL) %>% 
	select(subfamily,FL,Non_FL) %>% melt() %>% 
	mutate(subfamily_new = paste0(subfamily,"_",variable)) %>% 
	select(subfamily_new,value)

factor <- df_combine %>% group_by(cellLine) %>% summarize(sum = n()) %>% 
	mutate(factor = log2(sum/7818 + 1))

plotdf <- df_combine %>% 
	mutate(lengthGroup = ifelse(length >= 5999, "FL","Non_FL")) %>% 
	group_by(cellLine,V4,lengthGroup) %>% 
	summarise(count = n()) %>% 
	mutate(subfamily = factor(V4,levels = LINE1factor$V1)) %>% 
	mutate(subfamily_new = paste0(subfamily,"_",lengthGroup)) %>% 
	left_join(L1_cn_new) %>% 
	left_join(factor) %>% 
	mutate(norCount = count / factor) %>% 
	filter(value > 10) %>% 
	mutate(percent = norCount /value)


tmp <- plotdf %>% 
	filter(lengthGroup== "FL") %>% 
	dcast(cellLine ~ subfamily,value.var = "percent")
tmp[is.na(tmp)] <- 0
fwrite(tmp,"norm_FL_include_percent.tsv",sep = '\t' )
mat <- tmp[,-1]
rownames(mat) <- tmp$cellLine
col_fun = circlize::colorRamp2(c(0,0.03,0.07,0.09,0.1,0.13,1), colorRampPalette(rev(brewer.pal(n = 7, name =	"RdYlBu")))(7))
#col_fun = circlize::colorRamp2(c(0,0.05,0.1,0.2,1), colorRampPalette(rev(brewer.pal(n = 7, name =	"RdYlBu")))(5))
annotation_col <- data.frame(Group = rep(c("FL"),each=ncol(mat)))
label_L1 <- c("L1HS","L1PA2","L1PA3","L1PA4")
L1_position <- match(label_L1,colnames(mat))
rowanno = rowAnnotation(foo = anno_mark(at = L1_position, 
																				labels = label_L1,labels_gp = gpar(fontsize = 6)))
p1 <- ComplexHeatmap::pheatmap(t(mat),
															 col = col_fun,
															 border =F,cellheight = 1,cellwidth = 5.5,
															 right_annotation = rowanno,
															 annotation_row = annotation_col,
															 cluster_cols = F,fontsize = 5,show_rownames = F,
															 cluster_rows = F)
p1
tmp <- plotdf %>% 
	filter(lengthGroup== "Non_FL") %>% 
	dcast(cellLine ~ subfamily,value.var = "percent")
tmp[is.na(tmp)] <- 0
fwrite(tmp,"norm_nonFL_include_percent.tsv",sep = '\t' )
mat <- tmp[,-1]
rownames(mat) <- tmp$cellLine
col_fun = circlize::colorRamp2(c(0,0.1,0.2,0.3,1), colorRampPalette(rev(brewer.pal(n = 7, name =	"RdYlBu")))(5))
annotation_col <- data.frame(Group = rep(c("Non-FL"),each=ncol(mat)))
label_L1 <- c("L1PBa1","L1P2","L1PBa","L1P3b")
L1_position <- match(label_L1,colnames(mat))
rowanno = rowAnnotation(foo = anno_mark(at = L1_position, 
																				labels = label_L1,labels_gp = gpar(fontsize = 6)))
p2 <- ComplexHeatmap::pheatmap(t(mat),
															 annotation_row = annotation_col,
															 col = col_fun,cellheight = 1,cellwidth = 5.5,
															 border = F,fontsize = 5,
															 right_annotation = rowanno,
															 cluster_cols = F,show_rownames = F,
															 cluster_rows = F)
pdf("heatmap_percent_FullLength.pdf",width = 3,height = 5)
p <- p1%v%p2
draw(p,ht_gap = unit(1, "mm"))
dev.off()
