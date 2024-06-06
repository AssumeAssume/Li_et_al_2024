#!/analysis2/software/miniconda3/bin/Rscript
library(tidyverse)
library(magrittr)
library(data.table)
library(glue)
library(ggpubr)
theme_set(theme_bw())
rm(list=ls())
filter <- dplyr::filter
select <- dplyr::select
desc <- dplyr::desc
reduce <- purrr::reduce
rename <- dplyr::rename
readr::local_edition(1)
options(stringsAsFactors = FALSE)
set.seed(629)
source("/analysis/lixf/my_script/my_R_function.R")



setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee1/major_point6_CRISPRai_L1_proximalGene/NonFusion_K27ac_increase_Allk27ac_ContactGene")

distanceCutoff <- 32

DESeq2.knowngene.CTBP1_cKO <- fread("/analysis/lixf/CRISPRa_i/NCCIT/RNA/CRISPRi/5-DESeq2/KnownGene/Combine_add_KnownGene.Deseq2.L1_mix_vs_SAFE.result.txt") 

file <- dir("/analysis/lixf/CRISPRa_i/NCCIT/RNA/CRISPRi/7-windowBed/NonFusion_K27ac_increase_Allk27ac/",pattern = "stream_genename.txt",full.names = T)
file
df <- map(file,~fread(.x)) %>% reduce(.,rbind)

colnames(df) <- c("rowname","distance")
##### 
fusionGeneName <- fread("/analysis/lixf/CRISPRa_i/NCCIT/RNA/CRISPRi/get_fusion_L1/All_K27ac/Promoter_L1_fusion.gene.name.txt",header = F) 
Loop <- fread("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee1/major_point6_CRISPRai_L1_proximalGene/DeepLoop_intersect/K27ac_inorder_L1_Gene_pair.txt",header = F)
LoopGene <- Loop$V10 %>% unique

FC_df <- left_join(df,DESeq2.knowngene.CTBP1_cKO) %>% 
	mutate(GeneGroup  = ifelse(rowname %in% LoopGene,"Loop Gene","Non-loop gene"))

plotdf <- FC_df %>% drop_na() %>%
	filter(! rowname %in% fusionGeneName$V1)  


p2 <- plotdf %>% 
	ggplot(aes(x = GeneGroup, y =log2FoldChange))+
	geom_boxplot(outlier.alpha = 0,width = 0.5)+
	scale_y_continuous(expand = c(0,0),breaks = seq(-1,0.5,0.5))+
	coord_cartesian( ylim = c(-1,0.85))+
	geom_hline(yintercept = 0,color = "black",alpha = 0.7,linetype = 2)+
	stat_compare_means(label.y = c(0.48),method = "wilcox.test",size= 6/ggplot2::.pt)+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "right",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
	)+guides(fill=guide_legend(title="")) 
p2


# CRISPRa -----------------------------------------------------------------------------------------------------------------------------


DESeq2.knowngene.CTBP1_cKO <- fread("/analysis/lixf/CRISPRa_i/NCCIT/RNA/CRISPRa/5-DESeq2/KnownGene/Combine_add_KnownGene.Deseq2.L1_mix_vs_SAFE.result.txt") 

file <- dir("/analysis/lixf/CRISPRa_i/NCCIT/RNA/CRISPRa/7-windowBed/NonFusion_K27ac_increase_Allk27ac",pattern = "stream_genename.txt",full.names = T)
df <- map(file,~fread(.x)) %>% reduce(.,rbind)
colnames(df) <- c("rowname","distance")
##### 
fusionGeneName <- fread("/analysis/lixf/CRISPRa_i/NCCIT/RNA/CRISPRa/get_fusion_L1/All_K27ac/Promoter_L1_fusion.gene.name.txt",header = F) 
Loop <- fread("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee1/major_point6_CRISPRai_L1_proximalGene/DeepLoop_intersect/K27ac_inorder_L1_Gene_pair.txt",header = F)
LoopGene <- Loop$V10 %>% unique
FC_df <- left_join(df,DESeq2.knowngene.CTBP1_cKO) %>% 
	mutate(GeneGroup  = ifelse(rowname %in% LoopGene,"Loop Gene","Non-loop gene")) 


plotdf <- FC_df %>% drop_na() %>%
	filter(! rowname %in% fusionGeneName$V1)

plotdf %>% group_by(GeneGroup) %>% count()
pp2 <- plotdf %>% 
	ggplot(aes(x = GeneGroup, y =log2FoldChange))+
	geom_boxplot(outlier.alpha = 0,width = 0.5)+
	scale_y_continuous(expand = c(0,0),breaks = seq(-2,3,1))+
	coord_cartesian( ylim = c(-2.5,3.15))+
	geom_hline(yintercept = 0,color = "black",alpha = 0.7,linetype = 2)+
	stat_compare_means(label.y = c(1.38),method = "wilcox.test",size= 6/ggplot2::.pt)+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "right",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
	)+guides(fill=guide_legend(title="")) 
pp2

pCombine2 <- (pp2+p2)+plot_layout(guides = 'collect')&theme(panel.grid.major = element_blank(),
																														panel.grid.minor = element_blank(),
																														text=element_text(size= 8,family = "ArialMT",color = "black"),
																														plot.subtitle = element_text(hjust = 1),
																														axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
																														plot.margin = margin(0.3,1,0.3,0.3, "cm"),
																														strip.placement = "outside",
																														strip.background = element_blank(),
																														legend.position = "right",
																														legend.key.width= unit(6, 'pt'),
																														legend.key.height= unit(4, 'pt'),
																														axis.ticks = element_line(colour = "black", size = 0.25),
																														axis.ticks.length=unit(1.5, "pt")
)
pCombine2
ggsave(plot = pCombine2,glue("boxplot_combined.pdf"),width = 3.2,height = 2.1)
