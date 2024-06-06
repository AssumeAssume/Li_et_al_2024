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


setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure5_NCCIT_CRISPRai/CRISPRai_gene_expression_v2/")
CRIPSRa <- fread("/analysis/lixf/CRISPRa_i/NCCIT/RNA/CRISPRa/5-DESeq2/KnownGene/Combine_add_KnownGene.Deseq2.L1_mix_vs_SAFE.result.txt")
CRIPSRa[is.na(CRIPSRa)] <- 1
CRIPSRi <- fread("/analysis/lixf/CRISPRa_i/NCCIT/RNA/CRISPRi/5-DESeq2/KnownGene/Combine_add_KnownGene.Deseq2.L1_mix_vs_SAFE.result.txt")
CRIPSRi[is.na(CRIPSRi)] <- 1


tmp1 <- CRIPSRa %>% select(rowname,log2FoldChange,padj,pvalue ) %>% rename(`log2FC CRISPRa` = log2FoldChange,Apadj = padj,Apvalue = pvalue)
tmp2 <- CRIPSRi %>% select(rowname,log2FoldChange,padj,pvalue) %>% rename(`log2FC CRISPRi` = log2FoldChange,Ipadj = padj,Ipvalue = pvalue)
combine_df <- full_join(tmp1,tmp2)
combine_df$Apadj[is.na(combine_df$Apadj)] <- 1
combine_df$Ipadj[is.na(combine_df$Ipadj)] <- 1
combine_df$`log2FC CRISPRa`[is.na(combine_df$`log2FC CRISPRa`)] <- 0
combine_df$`log2FC CRISPRi`[is.na(combine_df$`log2FC CRISPRi`)] <- 0
###

combine_df <- combine_df


colorManual <- c("L1-regulated" = '#e41a1c',others = "grey")
sizeManual <- c("L1-regulated" = 0.4,others = 0.15)
alphaManual <- c("L1-regulated" = 0.9,others = 0.5)



tmpp <- fread("/analysis/lixf/CRISPRa_i/NCCIT/RNA/CRISPRa/closest_L1/padj0.1_lfc0_L1_regulated_gene_nearestL1_v2.txt",header = T)
tmpp <- tmpp %>% filter(`log2FC CRISPRa` > 0) 
colnames(tmpp)[1] <- "V1"
loosen_number <- nrow(tmpp)

combine_df %>% 
	mutate(group = ifelse(rowname %in% tmpp$V1, "L1-regulated","others")) %>% 
	ggplot(aes( x = `log2FC CRISPRa`, y = `log2FC CRISPRi`,color = group, size = group , alpha = group))+
	geom_point(stroke = 0.1)+
	geom_hline(yintercept = 0,color = "black",alpha = 0.7,linetype = 2,size = 0.3)+
	geom_vline(xintercept = 0,color = "black",alpha = 0.7,linetype = 2,size = 0.3)+
	scale_color_manual(values = colorManual,name = "")+
	scale_size_manual(values = sizeManual)+
	scale_alpha_manual(values = alphaManual)+
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
	)+
	annotate("text", x = 4.5, y = -1.75, label = glue("Number of loosen-criteria L1 regulated genes : \n {loosen_number}"),size = 6/ggplot2::.pt)+
	guides(color=guide_legend(title="family",override.aes = list(size=1.2)),alpha="none",size="none",text = "none")
ggsave("CRISPRai_gene_expression_added_genesets.pdf",height=2.5,width = 3.7,useDingbats = F)

