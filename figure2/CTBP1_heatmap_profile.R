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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/")
file <- system("tree -f -i /analysis/lixf/sgRNA/human/LINE_screening_combine//4-Deeptools/FLL1_profile/FL_L1PA_L1MA* |grep -i combine_profile_grouped.txt|grep -v 'FL_L1PA_L1MA_STARR/combine_profil'",intern = T)
color <- c("#bd3106", "#eebe04", "#c3d6ce" ,"#454b87")
plot_profile <- function(x,color = color){
	# x <- "/analysis/lixf/sgRNA/human/LINE_screening_co/4-Deeptools/FLL1_profile/FL_L1PA_L1MA_STARR/combine_profile_grouped.txt"
profile <- vroom::vroom(x,col_names=F,skip = 2)

colnames(profile)[1:2] <- c("sample","group")
meltd <- reshape2::melt(profile,id.vars=c("sample","group")) %>%  mutate(group = str_remove(group,".bed"))

meltd$group <- str_replace(meltd$group," \\(","\\\n\\(")
meltd$group <- fct_inorder(meltd$group)
### sample factor
#meltd$sample <- factor(meltd$sample,levels = c("CTBP1","PolII","K4me3","K9ac","K27ac","K9me3" ,"STARR"  ))

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
tmp <- meltd$sample %>% unique
if(length(tmp) == 2) meltd$sample <- factor(meltd$sample,levels = c(tmp[str_detect(tmp,"WT")],tmp[str_detect(tmp,"KO")]))
	
p <- ggplot(meltd, aes(x=variable, y=value, group=sample,color = sample)) + 
	geom_line(method = "loess",stat="smooth", se = FALSE,size=0.2, span=0.1,alpha = 0.8) +
	ylab(paste0("ChIP-Seq Signal (FPKM) ")) + 
	xlab("") +
	#geom_hline(yintercept = 0,color = "red",alpha = 0.7,linetype = 2,size = 0.2)+
	scale_x_discrete(breaks=pos, labels= xlabels) + 
	#scale_y_continuous(breaks=seq(0,0.12,0.02))+
	scale_color_manual(values = c("#454b87","#bd3106")) +
	facet_wrap(~group,nrow = 1) + 
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
#theme(strip.text = element_text(size = rel(0.5), color="black"), strip.background = element_rect(fill = 'white', colour = 'white'))
}
str_remove_all(file,".*_L1MA_|/combine.*")
plotList <- map(file,function(x){plot_profile(x)+ggtitle(str_remove_all(x,".*_L1MA_|/combine.*"))})
#pp <- reduce(plotList,`/`)
#file
file
pp <- plotList[[6]]/plotList[[1]]/plotList[[2]]/plotList[[3]]/(plotList[[5]]+scale_y_continuous(limits =c(0,0.05)))/plotList[[4]] &theme(legend.position = "right")# +
	#plot_layout(height = c())
pp <- reduce(plotList,`/`)&theme(legend.position = "right")

ggsave(plot = pp,"20221221_RNA_histone_marker_along_CTBP1_boundL1.pdf",width = 7.5,height = 12.5)
