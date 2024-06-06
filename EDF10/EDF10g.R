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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/RC_ASO_plot/")
ExpDf <- fread("/analysis/lixf/CRISPRa_i/mESC/RNA/ASO2c/5-DESeq2/KnownGene/Combine_Deseq2.RC_vs_SCR.result.txt")%>% mutate(sample = "L2C RC")
ExpDf2 <- fread("/analysis/lixf/CRISPRa_i/mESC/RNA/ASO2c/5-DESeq2/KnownGene/Combine_removeL1_1Deseq2.L1_vs_SCR.result.txt")%>% mutate(sample = "L2C ASO")

geneSets <- fread("/analysis/lixf/CRISPRa_i/mESC/RNA/public_gene_sets/2c_genesets.txt") %>% rename(symbol = GSE33923) %>%  mutate(geneset = "2C specific genes")  %>% select(geneset,symbol) %>% distinct()

p1 <- ExpDf %>% 
	mutate(group = ifelse(rowname %in% geneSets$symbol,"2C-specific","others")) %>% 
	filter(group == "2C-specific") %>% 
	arrange(desc(group)) %>% 
	ggplot(aes(x = log10(baseMean),y = log2FoldChange,color = group,size = group))+
	geom_point(stroke= 0.1)+
	scale_color_manual(values = c("2C-specific"='red',"others"='grey'))+
	scale_size_manual(values = c("2C-specific"=0.9,"others"=0.5))+
	geom_hline(yintercept = 0,color = "red",alpha = 0.7,linetype = 2)+
	ggtitle(" 2-cell, RC-L1 vs SCR")+
	ylim(-3,3)


p2 <- ExpDf2 %>% 
	mutate(group = ifelse(rowname %in% geneSets$symbol,"2C-specific","others")) %>% 
	arrange(desc(group)) %>% 
	filter(group == "2C-specific") %>% 
	ggplot(aes(x = log10(baseMean),y = log2FoldChange,color = group,size = group))+
	geom_point(stroke= 0.1)+
	scale_color_manual(values = c("2C-specific"='red',"others"='grey'))+
	scale_size_manual(values = c("2C-specific"=0.9,"others"=0.5))+
	geom_hline(yintercept = 0,color = "red",alpha = 0.7,linetype = 2)+
	ggtitle(" 2-cell, L1 vs SCR")+	ylim(-3,3)#+xlim(0,20)
(p1+p2)+plot_layout(guides='collect')&theme(panel.grid.major = element_blank(),
																						panel.grid.minor = element_blank(),
																						text=element_text(size= 8,family = "ArialMT",color = "black"),
																						plot.subtitle = element_text(hjust = 1),
																						axis.text.x = element_text(angle=0,vjust=0.5,hjust = 0),
																						plot.margin = margin(0.3,1,0.3,0.3, "cm"),
																						strip.placement = "outside",
																						strip.background = element_blank(),
																						legend.position = "bottom",
																						legend.key.width= unit(6, 'pt'),
																						legend.key.height= unit(4, 'pt'),
																						axis.ticks = element_line(colour = "black", size = 0.25),
																						axis.ticks.length=unit(1.5, "pt"),
																						axis.line =  element_line(colour = "black"),
																						
)
ggsave("MAplot_gene_exp_RC_L1AMO.pdf",width = 5.5,height = 3,useDingbats = F)
