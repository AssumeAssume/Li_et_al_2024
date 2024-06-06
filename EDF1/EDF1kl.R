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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/overlap_with_Other_Screen/overlap_with_Liu_GO")

df_filter <- fread("combine_TopGO.txt")

# supp --------------------------------------------------------------------------------------------------------------------------------


df_supp <- df_filter %>% filter(group2 == "suppressor")

df_term <- df_supp %>% 
	mutate(Term = str_remove(Term," \\(GO.*"),
				 Term=fct_reorder(Term, -log10(`P-value`),mean))
p1 <- df_term %>% 
	ggplot(aes(x = group, y = Term,color = -log10(`P-value`),size = overlap_count))+
	geom_point()+
	facet_wrap(~group2,scales = "free")+
	#scale_x_discrete(expand = c(0.01,0.1))+
	#scale_y_discrete(expand = c(0,0.5))+
	scale_color_gradient(high = "#2980B9",
											 low = "#aee3f9",
											 name= "-log10 P value")+
	scale_size_continuous(name = "Counts",breaks = c(5,10,20,30),range = c(0.3,3.8),limits = c(0,30))+
	theme_classic()+
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
	)+xlab("")+ylab("")


# acti --------------------------------------------------------------------------------------------------------------------------------


df_act <- df_filter %>% filter(group2 == "activator")

df_term <- df_act %>% 
	mutate(Term = str_remove(Term," \\(GO.*"),
				 Term=fct_reorder(Term, -log10(`P-value`),mean))
p2 <- df_term %>% 
	ggplot(aes(x = group, y = Term,color = -log10(`P-value`),size = overlap_count))+
	geom_point()+
	facet_wrap(~group2,scales = "free")+
	#scale_x_discrete(expand = c(0.01,0.1))+
	#scale_y_discrete(expand = c(0,0.5))+
	scale_color_gradient(high = "#2980B9",
											 low = "#aee3f9",
											 name= "-log10 P value")+
	scale_size_continuous(name = "Counts",breaks = c(5,10,20,30),range = c(0.3,3.8),limits = c(0,30))+
	theme_classic()+
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
	)+xlab("")+ylab("")
(p1/p2)+plot_layout(guides = "collect")
ggsave("combine_GOcompare_filter_count.pdf",width = 5.5,height = 3.5,useDingbats = F)


