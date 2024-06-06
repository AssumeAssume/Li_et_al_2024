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


setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure2/KMT2D/absolute_distance/")
x <- "./K27ac_increased_L1/distance.bed"

tmp <- fread(x) %>% 
	filter(V13 != -1)

tmp2 <-  table(cut_width(tmp$V13, 10000,boundary = 0,label =FALSE)) %>% as.data.frame() %>% mutate(group = "activated L1s in KMT2D KO")


# background --------------------------------------------------------------------------------------------------------------------------


file <- dir("./background_exp","distance",full.names = T)

background <- map(file, ~ fread(.x) %>% filter(V13 != -1)) %>% reduce(rbind)

background_bin <-  table(cut_width(background$V13, 10000,boundary = 0,label =FALSE)) %>% as.data.frame()
normalized <- background_bin %>% mutate(Freq = Freq/1000) %>% mutate(group = "background expectation")

color <- c("#7BB6A5","black")
df_combine <- rbind(tmp2,normalized) 
fwrite(df_combine,"./up_gene_distance_df.tsv",sep = '\t')
p <- df_combine %>% 
	mutate(Var1 =as.numeric(Var1)) %>% 
	filter(Var1 <= 100) %>% 
	mutate(Var1 = Var1 * 10) %>%
	ggplot(aes(x = Var1,y = Freq,group = group,color = group))+
	coord_cartesian(xlim = c(0,500))+
	geom_smooth(span = 0.1,se = F,size = 0.5)+
	scale_x_continuous(expand = c(0,0))+
	scale_color_manual(values = (color),name = "")+
	theme_classic()+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=0,vjust=0),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "bottom",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt")
	)+xlab("Distance (kb) from derepressed genes")+ylab("L1 elements")
p
ggsave(plot = p ,"k27ac_changedL1_plot.pdf",width = 1.8,height = 1.8)

