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
setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/Figure2g_add_down_gene")

down <- fread("./down_gene_distance_df.tsv")
up <- fread("up_gene_distance_df.tsv") %>% mutate(group = ifelse(str_detect(group,"K27ac"),"up-regulated genes","up background expectation"))

df_combine <- rbind(down,up) %>% mutate(linetype = ifelse(str_detect(group,"down"),"2","1"))
color <- c("down background expectation"="grey","down-regulated genes"="#7BB6A5",
					 "up background expectation"="black","up-regulated genes"="#B33317")

conf_1 <- fread("./down_background_Conf.tsv")
conf_2 <- fread("./up_background_Conf.tsv")

p <- df_combine %>% 
	mutate(Var1 =as.numeric(Var1)) %>% 
	mutate(Var1 = Var1 * 10) %>%
	ggplot(aes(x = Var1))+
	coord_cartesian(xlim = c(0,500))+
	geom_smooth(aes(y = Freq,group = group,color = group),span = 0.1,se = F,size = 0.2)+
	geom_smooth(data = conf_2 , aes(x = Var1*10, y = V1),span=0.1,se = F, alpha = 0.5,linetype="dotted",color="black",size = 0.2)+
	geom_smooth(data = conf_2 , aes(x = Var1*10, y = V2),span = 0.1,se = F, alpha = 0.5,linetype="dotted",color="black",size = 0.2)+
	geom_smooth(data = conf_1 , aes(x = Var1*10, y = V1),span=0.1,se = F, alpha = 0.5,linetype="dotted",color="grey",size = 0.2)+
	geom_smooth(data = conf_1 , aes(x = Var1*10, y = V2),span = 0.1,se = F, alpha = 0.5,linetype="dotted",color="grey",size = 0.2)+
	
	#geom_ribbon(data = conf_2 , aes(x = Var1*10, ymin = V1, ymax = V2), alpha = 0.7,linetype="dashed",color="black",size = 0.3)+
	#geom_ribbon(data = conf_1 , aes(x = Var1*10, ymin = V1, ymax = V2), alpha = 0.7,linetype="dashed",color="grey",size = 0.3)+
	#geom_point()+
	scale_x_continuous(expand = c(0,0))+
	scale_color_manual(values = rev(color),name = "")+
	theme_classic()+
	guides(
		linetype = "none",
		color=guide_legend(nrow=2,byrow=TRUE)
	)+
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
ggsave(plot = p ,"k27ac_changedL1_plot_with_down_add_95Line.pdf",width = 2.4,height = 1.8)

