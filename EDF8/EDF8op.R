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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/ZGA_pairwise_correlation/")

zga_l1_regulated <- fread("/analysis2/lixf/evolution/RNA/rhesus_embryo/Wang_GenomeResearch_2017_rhesus_embryo/5-DESeq2/TPM/zga_gene_define/n150_gene_name.txt",header = F)
L1_regulated_gene_added <- fread("/analysis/lixf/CRISPRa_i/NCCIT/RNA/CRISPRa/closest_L1/padj0.1_lfc0_L1_regulated_gene_nearestL1_v2.txt")
ZGA_L1_df <- L1_regulated_gene_added %>% 
	filter(rowname%in% zga_l1_regulated$V1) %>% 
	mutate(L1UniqueID= paste(L1Subfamily,L1Chr,L1Start+1,L1End,L1Strand,sep = '_'))
CRISPRa_L1 <- fread("/analysis/lixf/CRISPRa_i/NCCIT/RNA/CRISPRa/5-DESeq2/RepeatMasker_Unique/Combine_add_RepeatMasker_Unique.Deseq2.L1_mix_vs_SAFE.result.txt") %>% 
	filter(rowname %in% ZGA_L1_df$L1UniqueID)
colnames(CRISPRa_L1)[1] <- "L1UniqueID"

histone <- fread("/analysis/lixf/CRISPRa_i/NCCIT/ChIP/histone_enrichment/histon_enrichment_individual.txt")
promoter <- histone %>% 
	mutate(L1UniqueID = paste(L1_subfamily,chr,start,end,strand,sep ='_')) %>%
	#filter(Enrich_H3K4me3_WT != 0 & Enrich_PolII_WT !=0 & Enrich_H3K4me3_KO!=0 & Enrich_PolII_KO != 0) %>% 
	mutate(LFCPolII_A= log2((Enrich_PolII_CRISPRa+0.01)/(Enrich_PolII_CRISPRa_NC+0.01)),
				 LFCPolII_I= log2((Enrich_PolII_CRISPRi+0.01)/(Enrich_PolII_CRISPRi_NC+0.01)),
				 LFCK27ac_A= log2((Enrich_H3K27ac_CRISPRa+0.01)/(Enrich_H3K27ac_CRISPRa_NC+0.01)),
				 LFCK27ac_I= log2((Enrich_H3K27ac_CRISPRi+0.01)/(Enrich_H3K27ac_CRISPRi_NC+0.01)),
				 LFCCTCF_A= log2((Enrich_CTCF_CRISPRa+0.01)/(Enrich_CTCF_CRISPRa_NC+0.01)),
				 LFCRAD21_A= log2((Enrich_RAD21_CRISPRa+0.01)/(Enrich_RAD21_CRISPRa_NC+0.01)),
	) %>% drop_na() %>% 
	select(L1UniqueID,contains("LFC"))



pairwise <- left_join(ZGA_L1_df,promoter) %>% drop_na()

pairwiseA <- pairwise

model <- lm(`log2FC CRISPRa`~LFCK27ac_A,pairwiseA )
modelSummary <- glance(model) %>% 
	mutate(R = sqrt(r.squared))
modelSummary
modelSummary$R
N = pairwiseA %>% nrow()

p <- pairwiseA %>% 
	ggplot(aes(x = LFCK27ac_A,y = `log2FC CRISPRa`))+
	geom_point(alpha = 0.9,size = 0.9,stroke = 0.1)+
	geom_smooth(method = "lm",size = 0.5)+
	annotate(geom = "text", x = 1.5,y = 3.0,label = glue("r = {round(modelSummary$R,3)}\n P = {signif(modelSummary$p.value,digits =2)}\n N={N} "),
					 size = 6/ggplot2::.pt)+
	xlab("L1 K27ac log2[CRISPRa/Ctrl] ")+ylab("ZGA Gene expression \n log2[CRISPRa/Ctrl]")+
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
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt")
	)



pairwiseI <- pairwise 
model <- lm(`log2FC CRISPRi`~LFCK27ac_I,pairwiseI)
modelSummary <- glance(model) %>% 
	mutate(R = sqrt(r.squared))
modelSummary
modelSummary$R
N = pairwiseI %>% nrow()

p2 <- pairwiseI %>% 
	ggplot(aes(x = LFCK27ac_I,y = `log2FC CRISPRi`))+
	geom_point(alpha = 0.9,size = 0.9,stroke = 0.1)+
	geom_smooth(method = "lm",size = 0.5)+
	annotate(geom = "text", x = -0.5,y = -0.3,label = glue("r = {round(modelSummary$R,3)}\n P = {signif(modelSummary$p.value,digits =2)}\n N={N} "),
					 size = 6/ggplot2::.pt)+
	xlab("L1 K27ac log2[CRISPRi/Ctrl] ")+ylab("ZGA Gene expression \n log2[CRISPRi/Ctrl]")+
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
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt")
	)


p+p2
ggsave("combine_CRISPRai_ZGA_pairwise_K27ac.pdf",width = 7,height = 2)

