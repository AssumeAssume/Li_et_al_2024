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
setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/intron_splicing/intron_count/")
intron <- fread("/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/intron_splicing/intron_count/multiqc_data_3///multiqc_rseqc_read_distribution.txt")
intronDf <- intron %>% 
	select(Sample,total_tags,introns_tag_count,introns_tag_pct) %>% 
	mutate(intron_fraction = introns_tag_count/total_tags) %>% 
	mutate(group = ifelse(str_detect(Sample,"NC|Control"),"Control",ifelse(str_detect(Sample,"ataxia"),"Ataxia","KO")),
				 group2 = ifelse(str_detect(Sample,"ataxia|Control"),"Takahashi et al.","this paper")
	) %>% 
	mutate(gene = str_remove_all(Sample,"_shRNA[0-9]|_CRISPR.*$|NC_|_rep[0-9]|-[0-9]{1,2}|_[0-9]|[0-9]{3}-|clonal_"),
				 gene = ifelse(str_detect(Sample,"ataxia|Control"),"Ataxia",gene)
	) %>% 
	separate_rows( gene, sep = '_', convert = TRUE)
intronDf %>% 
	mutate(gene = fct_reorder(gene,introns_tag_pct,mean)  ) %>% 
	mutate(group2 = factor(group2,c("this paper","ENCODE","Takahashi et al."))) %>% 
	#filter(group2 == "this paper") %>% 
	ggplot(aes(x = gene,y = introns_tag_pct/100,color = group))+
	geom_point(alpha = 0.8,size = 1,stroke = 0.1)+
	scale_y_continuous(limits = c(0,1),label = percent,breaks = c(0,0.2,0.4,0.6,0.8,1))+
	ylab("intronic read percentage")+
	scale_color_jco()+
	facet_grid(~group2,space="free",scales = "free")+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "bottom",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt"),
				axis.line =  element_line(colour = "black", size = 0.25)
	)
ggsave("intron_pct_with_ATAXIA_color2.pdf",width = 4,height =2.5)

intronDf %>% 
	mutate(gene = fct_reorder(gene,introns_tag_pct,mean)  ) %>% 
	mutate(group2 = factor(group2,c("this paper","ENCODE","Takahashi et al."))) %>% 
	filter(group2 == "this paper") %>% 
	ggplot(aes(x = gene,y = introns_tag_pct/100,color = group))+
	geom_point(alpha = 0.8,size = 1,stroke = 0.1)+
	scale_y_continuous(limits = c(0,0.5),label = percent,breaks = c(0,0.16,0.2,0.4,0.6,0.8,1))+
	ylab("intronic read percentage")+
	scale_color_jco()+
	geom_hline(yintercept = 0.16,color = "red",alpha = 0.7,linetype = 2)+
	facet_grid(~group2,space="free",scales = "free")+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "bottom",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt"),
				axis.line =  element_line(colour = "black", size = 0.25)
	)
ggsave("intron_pct_without_ATAXIA_color2.pdf",width = 4,height =2.5)
