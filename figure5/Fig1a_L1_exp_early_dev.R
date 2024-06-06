library(tidyverse)
library(magrittr)
library(data.table)
library(glue)
library(ungeviz)
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
setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/mouseCRISPRai/L1_exp_early_dev")

df <- fread("/analysis/public_data_analysis/mouse_embryo/Deng_2014_Science_GSE45719/7-ExpressionPlot/RepeatMasker_Unique.TPM_noFilter.tsv")
i="FL L1"
a <- df %>% 
	filter(str_detect(rowname,"^L1")) %>% 
	mutate(end = str_split_fixed(str_split_fixed(rowname,"_chr",2)[,2],"_",4)[,3] %>% as.numeric(),
				 start = str_split_fixed(str_split_fixed(rowname,"_chr",2)[,2],"_",4)[,2]%>% as.numeric()  )%>%
	filter(end - start >= 5000) %>% 
	select(-end,-start) 
N <- nrow(a)
b <- data.frame(value=apply(a[,-1],2,sum)) %>% rownames_to_column() %>% 
	mutate(rowname = str_replace(rowname,"16cell","morula"))
b$group <- factor( b$rowname %>% str_remove("_.*|[0-9]$"),levels=c("zy","early2cell","mid2cell" ,"late2cell" ,"4cell" ,"8cell","morula","earlyblast","midblast","lateblast"))
colorManual <- c("#5E251B","#A53222","#C78C3A","#6B4887","#5765A3","#AECDE2","#F0E759","#DAD95F","#5B7944")
# my_comparisons <- list(c("mid2cell","zy"),c("mid2cell","early2cell"),
# 											 c("late2cell","zy"),c("late2cell","early2cell"),
# 											 c("4cell","zy"),c("4cell","early2cell"),
# 											 c("mid2cell","8cell"),c("mid2cell","morula"),c("mid2cell","lateblast"),
# 											 c("late2cell","8cell"),c("late2cell","morula"),c("late2cell","lateblast"),
# 											 c("4cell","8cell"),c("4cell","morula"),c("4cell","lateblast")
# 											 )
wilcox.test(b %>% filter(group == "early2cell") %>% pull(value),b %>% filter(str_detect(group,"zy")) %>% pull(value),alternative = "greater")
wilcox.test(b %>% filter(group == "mid2cell") %>% pull(value),b %>% filter(str_detect(group,"zy|early2cell")) %>% pull(value))
wilcox.test(b %>% filter(group == "late2cell") %>% pull(value),b %>% filter(str_detect(group,"zy|early2cell")) %>% pull(value))
wilcox.test(b %>% filter(group == "4cell") %>% pull(value),b %>% filter(str_detect(group,"zy|early2cell")) %>% pull(value))
wilcox.test(b %>% filter(group == "mid2cell") %>% pull(value),b %>% filter(str_detect(group,"8cell|lateblast")) %>% pull(value))
wilcox.test(b %>% filter(group == "late2cell") %>% pull(value),b %>% filter(str_detect(group,"8cell|lateblast")) %>% pull(value))
wilcox.test(b %>% filter(group == "4cell") %>% pull(value),b %>% filter(str_detect(group,"8cell|lateblast")) %>% pull(value))
wilcox.test(b %>% filter(group == "early2cell") %>% pull(value),b %>% filter(str_detect(group,"8cell|lateblast")) %>% pull(value))

b %>%
	filter(!str_detect(group,"earlyblast|midblast")) %>% group_by(group) %>% count()

pFL <- b %>%
	filter(!str_detect(group,"earlyblast|midblast")) %>% 
	ggplot(aes(x=group,y=value,color = group))+
	#geom_boxplot(aes(fill=group))+
	ggbeeswarm::geom_quasirandom(size = 0.5,stroke = 0.1)+
	stat_summary(fun = "mean", colour = "black",width = 0.5,size=0.3,geom = "hpline",alpha=0.8)+
	#facet_wrap(~group,nrow = 1,scales = "free_x")+
	scale_color_manual(values = colorManual)+
	#ggpubr::stat_compare_means(comparisons = my_comparisons,label =  'p.signif')+
	theme_classic()+
	scale_y_continuous(breaks = seq(0,100,25),limits = c(0,100))+
	labs(x="",y="TPM",title=paste0(i," N =(",N, ")"))+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "none",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt"),
				axis.line =  element_line(colour = "black", size = 0.25)
	)
	

pFL_L2 <- b %>%
	filter(!str_detect(group,"earlyblast|midblast")) %>% 
	ggplot(aes(x=group,y=log2(value),color = group))+
	#geom_boxplot(aes(fill=group))+
	ggbeeswarm::geom_quasirandom(size = 0.5,stroke = 0.1)+
	stat_summary(fun = "mean", colour = "black",width = 0.5,size=0.3,geom = "hpline",alpha=0.8)+
	#facet_wrap(~group,nrow = 1,scales = "free_x")+
	scale_color_manual(values = colorManual)+
	#ggpubr::stat_compare_means(comparisons = my_comparisons,label =  'p.signif')+
	theme_classic()+
	#scale_y_continuous(breaks = seq(0,100,25),limits = c(0,100))+
	labs(x="",y="Log2(TPM)",title=paste0(i," N =(",N, ")"))+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "none",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt"),
				axis.line =  element_line(colour = "black", size = 0.25)
	)
pFL
ggsave(plot = pFL+pFL_L2,"FL_L1_exp_with_P_noFilter_combined.pdf",width = 5,height = 2.5,useDingbats =F)
ggsave(plot = pFL,"noLog2_FL_L1_exp_with_P_noFilter_combined.pdf",width = 3,height = 2.5,useDingbats =F)
ggsave(plot = pFL_L2,"Log2_FL_L1_exp_with_P_noFilter_combined.pdf",width = 3,height = 2.5,useDingbats =F)
