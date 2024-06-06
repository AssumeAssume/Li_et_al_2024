library(tidyverse)
library(magrittr)
theme_set(theme_bw())
rm(list=ls())
filter <- dplyr::filter
select <- dplyr::select
desc <- dplyr::desc
reduce <- purrr::reduce
options(stringsAsFactors = FALSE)
set.seed(629)

setwd("/analysis2/lixf/DIS3/PAR_CLIP/4-Deeptools/DIS3_OVERLAP/FLL1_mappablity_combine/")
file <- c("/analysis2/lixf/DIS3/PAR_CLIP/4-Deeptools/DIS3_OVERLAP/FLL1_mappablity_combine/DIS3_profile_grouped.txt")

color <- c("#bd3106", "#eebe04", "#c3d6ce" ,"#454b87")
x <- file

profile <- vroom::vroom(x,col_names=F,skip = 2)

colnames(profile)[1:2] <- c("sample","group")
meltd <- reshape2::melt(profile,id.vars=c("sample","group")) %>%  mutate(group = str_remove(group,".bed"))

meltd$group <- str_replace(meltd$group," \\(","\\\n\\(")
meltd$group <- fct_inorder(meltd$group)

meltd$sample <- factor(meltd$sample,levels = c("PAR_DIS3","mappability"  ))

#### xlabels and label position
xlabels <- c("-1","L1 5'","L1 3'","1Kb")
xaxis_level <- levels(meltd$variable)
xaxis_length <- length(xaxis_level)
pos <-  c(xaxis_level[1], xaxis_level[xaxis_length*1/9],xaxis_level[xaxis_length*8/9],xaxis_level[xaxis_length])


p <- ggplot(meltd, aes(x=variable, y=value, group=sample)) + 
	geom_line(method = "loess",stat="smooth", se = FALSE,size=0.2, span=0.01,alpha = 0.8,color = "#454b87") +
	ylab(paste0("PAR-CLIP Signal (RPKM) ")) + 
	xlab("") +
	#geom_hline(yintercept = 0,color = "red",alpha = 0.7,linetype = 2,size = 0.2)+
	scale_x_discrete(breaks=pos, labels= xlabels) + 
	#scale_y_continuous(breaks=seq(0,0.12,0.02))+
	#scale_color_manual(values = c("#454b87")) +
	facet_wrap(~sample,ncol = 1,scales = "free_y") + 
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
				axis.ticks.length=unit(1.5, "pt"),
	)+guides(color=guide_legend(title="",nrow=1,byrow = T))+
	theme(panel.spacing.x=unit(0.1, "lines") )  ##### Alter just horizontal spacing between facets

ggsave(plot = p,"PARCLIP_profile.pdf",width = 2.7,height = 2.5)

p <- meltd %>% filter(sample == "mappability") %>% 
	ggplot(aes(x=variable, y=value, group=sample)) + 
	geom_line(method = "loess",stat="smooth", se = FALSE,size=0.2, span=0.01,alpha = 0.8,color = "#454b87") +
	ylab(paste0("PAR-CLIP Signal (RPKM) ")) + 
	xlab("") +
	#geom_hline(yintercept = 0,color = "red",alpha = 0.7,linetype = 2,size = 0.2)+
	scale_x_discrete(breaks=pos, labels= xlabels) + 
	#scale_y_continuous(breaks=seq(0,0.12,0.02))+
	#scale_color_manual(values = c("#454b87")) +
	#facet_wrap(~sample,ncol = 1,scales = "free_y") + 
	scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))+
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
				axis.ticks.length=unit(1.5, "pt"),
	)+guides(color=guide_legend(title="",nrow=1,byrow = T))+
	theme(panel.spacing.x=unit(0.1, "lines") )  ##### Alter just horizontal spacing between facets
p
ggsave("mappability.pdf",width = 2.7,height = 1.5)
