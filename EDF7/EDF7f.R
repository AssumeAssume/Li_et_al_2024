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
options(scipen = 200)
source("/analysis/lixf/my_script/my_R_function.R")


setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure5_NCCIT_CRISPRai")
exp <- fread("gene_expression_change.tsv")
distance <- fread("/analysis/lixf/CRISPRa_i/NCCIT/RNA/CRISPRa/closest_L1/assembled_gene_closest_minimum_selected_distance.txt",header = F)
TMP <- distance %>% select(V5,V13) %>% filter(!V13 == -1)
colnames(TMP) <- c("rowname","distance")
TMP <- TMP %>% 
	group_by(rowname) %>% 
	summarize(min_distance = min(distance)) %>% 
	ungroup() %>% 
	mutate(min_distance = min_distance/1000)
tmp <- left_join(exp,TMP) 
sum(is.na(TMP$min_distance))

combine_df <- left_join(exp,TMP)
nodistance_genenama <- combine_df[is.na(combine_df$min_distance)] %>% pull(rowname)
combine_df <- combine_df %>% filter(! rowname%in% nodistance_genenama) %>% 
	mutate(distance_bin = cut(min_distance,
														breaks = quantile(min_distance,probs = seq(0, 1, by = 0.1)), 
														labels = seq(1,10,1),
														include.lowest = TRUE
	),
	) 

combine_df %>% group_by(distance_bin) %>% count()
L1_regulated_gene_added <- fread("/analysis/lixf/CRISPRa_i/NCCIT/RNA/CRISPRa/closest_L1/padj0.1_lfc0_L1_regulated_gene_nearestL1_v2.txt")
geneName_lfc0 <- L1_regulated_gene_added %>% filter(`log2FC CRISPRa` > 0) %>% pull(rowname)
combine_df %>% 
	filter(rowname %in% geneName_lfc0) %>% 
	group_by(distance_bin) %>% 
	count()

tmpp <- combine_df$min_distance
levels_label = levels(cut(tmpp,breaks = quantile(tmpp,probs = seq(0, 1, by = 0.1),include.lowest = TRUE,right = F,dig.lab = 1)))

color <- colorRampPalette(c("orange", "blue"))(15)
color <- color[c(-2,-4,-6,-7,-8)] 
length(color)
combine_df$distance_bin %>% unique %>% sort
combine_df %>% 
	mutate(alpha = ifelse(as.numeric(distance_bin) < 2,  "whole", "alpha")) %>% 
	mutate(distance_bin = as.character(distance_bin)) %>% 
	arrange((distance_bin)) %>% 
	ggplot(aes( x = `log2FC CRISPRa`, y = `log2FC CRISPRi`,color = distance_bin))+
	geom_point(aes(alpha = alpha),size = 0.5,stroke = 0.2)+
	scale_alpha_manual(values = c("whole" = 1, "alpha" = 0.7))+
	geom_hline(yintercept = 0,color = "black",alpha = 0.7,linetype = 2,size = 0.3)+
	geom_vline(xintercept = 0,color = "black",alpha = 0.7,linetype = 2,size = 0.3)+
	#scale_color_gradient(low = "red",high = "blue")+
	scale_color_manual(values = color,name = "distance bin",	 guide = guide_coloursteps(show.limits = TRUE),
										 breaks = seq(1, 10, 1))+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=0,vjust=0),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "right",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
	)

ggsave("CRISPRai_gene_expression_padj0.1_geneDistance_Noalpha.pdf",height=2.5,width = 3.7,useDingbats = F)

