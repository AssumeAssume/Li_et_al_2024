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
setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/intron_splicing/SPLICEq/SE_stat/")
file <- dir("./",".tsv")
names(file) <- str_remove(file,".stat.tsv")
SE <- imap_dfr(file,~fread(.x),.id="Sample")
SE
encode <- fread("../../intron_count/ENCODE_redoBAM.txt",header = F) %>% 
	mutate(gene = str_remove_all(V1,"_shRNA[0-9]|_CRISPR.*$|NC_|_rep[0-9]|-[0-9]{1,2}|_[0-9]|[0-9]{3}-|clonal_"))
encode
Df <- SE %>% 
	mutate(group = ifelse(str_detect(Sample,"NC"),"Control","KO")) %>% 
	mutate(gene = str_remove_all(Sample,"_shRNA[0-9]|_CRISPR.*$|NC_|_rep[0-9]|-[0-9]{1,2}|_[0-9]|[0-9]{3}-|clonal_"),
				 group2= "this paper",
				 #group2 = ifelse(gene %in% encode$gene,"ENCODE","this paper")
	) %>% 
	separate_rows( gene, sep = '_', convert = TRUE)

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/intron_splicing/SPLICEq/SE_stat_ataxia/")
file <- dir("./",".tsv")
names(file) <- str_remove(file,".stat.tsv")
SE2 <- imap_dfr(file,~fread(.x),.id="Sample")
tmp <- SE2 %>% mutate(gene = "Ataxia",
											group = ifelse(str_detect(Sample,"Control"),"Control","Ataxia"),
											group2 = "Takahashi et al."
)
line <- tmp %>% filter(group == "NC") %>% group_by(group) %>% 
	summarize(mean=median(`mean(score)`))
Df2 <- rbind(Df,tmp)
Df2 %>% 
	mutate(gene = fct_reorder(gene,`mean(score)`,mean) %>% fct_rev ) %>% 
	mutate(group2 = factor(group2,c("this paper","ENCODE","Takahashi et al."))) %>% 
	arrange(group) %>% 
	ggplot(aes(x = gene,y = `mean(score)`,color = group))+
	geom_jitter(alpha = 0.8,size = 1,stroke = 0.1,width = 0.1)+
	geom_hline(yintercept = 0.96,color = "red",alpha = 0.7,linetype = 2)+
	scale_y_continuous(limits = c(0.8,1.001),label = percent)+
	facet_grid(~group2,space="free",scales='free')+
	ylab("splicing efficiency (%)")+
	scale_color_jco()+coord_cartesian(clip = T)+
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
ggsave("SE_rate_with_ataxia_color2.pdf",width = 4,height =2.5,useDingbats = F)

Df2 %>% 
	mutate(gene = fct_reorder(gene,`mean(score)`,mean) %>% fct_rev ) %>% 
	mutate(group2 = factor(group2,c("this paper","ENCODE","Takahashi et al."))) %>% 
	arrange(group) %>% 
	filter(group2 =="this paper") %>% 
	ggplot(aes(x = gene,y = `mean(score)`,color = group))+
	geom_jitter(alpha = 0.8,size = 1,stroke = 0.1,width = 0.1)+
	geom_hline(yintercept = 0.96,color = "red",alpha = 0.7,linetype = 2)+
	scale_y_continuous(limits = c(0.85,1.001),label = percent,breaks =c(0.85,0.9,0.95,1))+
	facet_grid(~group2,space="free",scales='free')+
	ylab("splicing efficiency (%)")+
	scale_color_jco()+coord_cartesian(clip = T)+
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
ggsave("SE_rate_without_ataxia_color2.pdf",width = 4,height =2.5,useDingbats = F)
