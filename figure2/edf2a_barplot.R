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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/second_rebuttle/Unique_Multiple_percent")
LINE1_factor <- fread("/LiuLab/reference/Human/GRCh38/TE/LINE_1_factor.txt",header = F) %>% unlist
unique <- fread("unique_count.rawcounts.tsv")
multiple <- fread("multiple_count.rawcounts.tsv")
df <- left_join(unique,multiple,by='Geneid')
colnames(df) <- c("subfamily","unique","all")
plotdf <- df	%>% mutate(multiple = all-unique) %>% 
	select(-all) %>% 
	melt() %>% 
	group_by(subfamily) %>% 
	mutate(percent = value/sum(value)) %>% 
	mutate(subfamily = factor(subfamily,LINE1_factor)
				 ) 
library(MetBrewer)
MetBrewer::colorblind_palettes
met.brewer("VanGogh3",7)
plotdf %>% 
	mutate(variable =fct_inorder(variable) %>% fct_rev) %>% 
	ggplot(aes(x = subfamily , y = percent,fill=variable))+
	geom_col()+
	scale_y_continuous(label = percent,expand = c(0,0.01),breaks = seq(0,1,0.2))+
	xlab("")+ylab("percent")+
	scale_fill_manual(values = c("#1A1919","#EE0000FF"),name="")+
	#scale_fill_aaas()
	#scale_fill_manual(values = met.brewer("VanGogh3",7)[c(1,3)],name="") +
	theme_classic()+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 6,family = "ArialMT",color = "black"),
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
	)+ylab("% of reads mapped to L1")
ggsave("reads_L1_percent.pdf",width = 3,height = 2)

