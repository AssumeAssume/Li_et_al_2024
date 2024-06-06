library(tidyverse)
library(magrittr)
library(patchwork)
library(foreach)
library(ungeviz)
theme_set(theme_bw())
rm(list=ls())
filter <- dplyr::filter
select <- dplyr::select
reduce <- purrr::reduce
options(stringsAsFactors = FALSE)
set.seed(629)


setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure2/DBR1//")


combine_df <- fread("./lariat_anno_filter_df.tsv")
combine_df$group <- factor(combine_df$group,levels=c("NC","DBR1 KD"))
p1 <- combine_df %>% 
		filter(X11 == "L1") %>% 
	  filter(group2 == "rRNA -") %>% 
		select(X1:X7,group2,group,unmapped,X4) %>% 
		distinct() %>% 
		group_by(group2,group,unmapped,X4) %>% 
		summarise(count = n()) %>% 
		mutate(CPMillion = count*1000000/unmapped) %>% 
		#filter(str_detect(group2,"DBR1")) %>% 
		ggplot(aes(x=group,y=CPMillion,color=group))+
		geom_point(size = 0.5)+
   	stat_summary(fun = "mean", colour = "black",width = 0.5,size=0.6,geom = "hpline",alpha=0.8)+
		ggsci::scale_color_aaas()+
		scale_y_continuous(limits = c(0,8))+
		#facet_wrap(~group2,scales = "free_x")+
	  ylab(" lariat counts / million \n unmapped singlets")+xlab("")+
	  ggtitle("lariats containing L1")+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "right",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt")
	)

# 

# total lariats -----------------------------------------------------------------------------------------------------------------------

normalize <- fread("./normalized_lariat_count.tsv")

normalize <- normalize %>% 
	mutate(group =  str_remove_all(group,"clonal|population|KD|_")) %>% 
	mutate(group = str_replace(group,"DBR1","DBR1 KD"),
				 group = factor(group,levels = c("NC","DBR1 KD"))) %>% 
	mutate(unmap_normalize_count = count*1000000 / (total_reads-mapped_reads))

p2 <- normalize %>% 
	filter(group2 == "rRNA -") %>% 
	ggplot(aes(x=group,y=unmap_normalize_count,color=group))+
	geom_point(size = 0.5)+
	stat_summary(fun = "mean", colour = "black",width = 0.5,size=0.6,geom = "hpline",alpha=0.8)+
	ggsci::scale_color_aaas()+
	#scale_y_continuous(limits = c(0,8))+
	#facet_wrap(~group2,scales = "free_x")+
	ylab(" lariat counts / million \n unmapped singlets")+xlab("")+
	ggtitle("All lariats")+
	scale_y_continuous(limits = c(0,130),breaks = seq(0,120,20))+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "right",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt")
	)
p2+p1+plot_layout(guides = "collect")
ggsave("lariats_combine.pdf",height=1.7,width = 3.5,useDingbats = F)
