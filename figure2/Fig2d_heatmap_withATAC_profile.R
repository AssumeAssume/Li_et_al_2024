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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure1_STARR_related/4-Deeptools/FLL1_profile/FLL1_with_ATAC_randomAssign_forProfile/")
color <- c("#bd3106", "#eebe04", "#c3d6ce" ,"#454b87")

x <- "./combine_profile_grouped.txt"
profile <- vroom::vroom(x,col_names=F,skip = 2)

colnames(profile)[1:2] <- c("sample","group")
meltd <- reshape2::melt(profile,id.vars=c("sample","group")) %>%  mutate(group = str_remove(group,".bed"))
levels <- meltd$sample %>% unique
meltd$group <- str_replace(meltd$group," \\(","\\\n\\(")
meltd$group <- fct_inorder(meltd$group)
### sample factor
meltd$sample <- factor(meltd$sample,levels =levels)

#### xlabels and label position
xlabels <- c("-3","L1 5'","L1 3'","3Kb")
xaxis_level <- levels(meltd$variable)
xaxis_length <- length(xaxis_level)
pos <-  c(xaxis_level[1], xaxis_level[xaxis_length/4],xaxis_level[xaxis_length*3/4],xaxis_level[xaxis_length])
#### color list
# col <- c("#464944", ###	WT hela
# 	"#5d605b","#babbb9", ### NC1
# 	"#0c1539","#2a3252","#49506b","#676d83","#868a9c", ### NC2
# 	"#e7b290","#ed9939","#b45d30","#f359da","#eb352a",### C1
# 	"#BFD8EC","#94A5D6","#2D73B6","#51B0E6","#033BF0" ### C6
# ) 
#### color list
#col <- color[1]
p <- ggplot(meltd, aes(x=variable, y=value, group = sample)) + 
	geom_smooth(method = "loess",stat="smooth", se = FALSE,size=0.1, span=0.1,alpha = 0.8) +
	ylab(paste0("Control Normalized Signal  ")) + 
	xlab("") +
	scale_x_discrete(breaks=pos, labels= xlabels)+
	#scale_y_continuous(breaks =seq(0.6,60,0.2))+
	facet_wrap(~sample,nrow = 1,scales = "free_y") + 
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

ggsave(plot = p,"STARRSeq_marker_along_FLL1.pdf",width = 7,height = 1.3)


p <- ggplot(meltd, aes(x=variable, y=value, group = sample)) + 
	geom_smooth(method = "loess",stat="smooth", se = FALSE,size=0.1, span=0.1,alpha = 0.8) +
	ylab(paste0("Control Normalized Signal  ")) + 
	xlab("") +
	scale_x_discrete(breaks=pos, labels= xlabels)+
	scale_y_continuous(breaks =seq(0.7,60,0.2))+
	facet_wrap(~sample,nrow = 1,scales = "free_y") + 
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

ggsave(plot = p,"STARRSeq_marker_along_FLL1_for_axis.pdf",width = 7,height = 9.3)




# FPKM --------------------------------------------------------------------------------------------------------------------------------

