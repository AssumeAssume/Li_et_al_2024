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

setwd("/analysis/lixf/tracks/STARR/mouse/Peng_GenomeBiology_2020/11-annoEnrichment_family_subfamily/subfamily_enrichment_5KL1_5UTR//")

file <- dir("./","ann.txt")
file
names(file) <- str_remove_all(file,"_enrichment.*|exp|_TF|_cell") 
annlist <- imap(file,~ fread(.x,skip = 6) %>% mutate(sample = .y)) %>% reduce(rbind)

####### filter out copies < 10, filter Promoter intron ..., only remain repeats, plot log2 Ratio (obs/exp)
family <- fread("/LiuLab/reference/Mouse/GRCm38/TE/TE_class/all.name",header = F)
colnames(family) <- c("subfamily","class","family")

# L1 ----------------------------------------------------------------------------------------------------------------------------------
LINE1factor <- fread("/LiuLab/reference/Mouse/GRCm38/TE/LINE_1_factor.txt",header = F) 

# remove promoter ---------------------------------------------------------------------------------------------------------------------

x <- annlist %>% 
	mutate(subfamily = str_remove(Annotation,"5UTR_FL_|non_FL_")) %>% 
	left_join(family) %>% 
	filter( class %in% c("LINE","SINE","DNA","LTR","Retroposon","Simple_repeat","Low_complexity")) %>% 
	filter( !str_detect(Annotation ,"\\?") ) %>% 
	filter(`Number of peaks` > 5 &  
				 	! Annotation == "ALR/Alpha" & 
				 	!str_detect(family,"\\?") 
	) 
labelLTR <- x %>% group_by(sample) %>% filter(class == "LTR") %>% slice_max(`Number of peaks`, n = 3) %>% ungroup
#label2 <- x %>% slice_min(`LogP enrichment (+values depleted)`, n =5)
#label_family <- x %>% group_by(class) %>% top_bottom(2,`Log2 Ratio (obs/exp)`,with_ties = F ) %>% ungroup
topL1 <- x %>% filter(str_detect(Annotation, "FL_")) %>% slice_min(`LogP enrichment (+values depleted)`, n = 20) %>%   pull(Annotation) %>% unique()
label_L1 <- x %>% filter(Annotation %in% topL1) 
#label_combine <- reduce(list(label_family,label,label2),rbind) %>% distinct()
label_combine <- reduce(list(label_L1),rbind) %>% distinct()
#colorManual <- c(LINE = '#A73030FF',SINE='#EFC000FF',LTR='#0073C2FF',Retroposon='#868686FF',"DNA"='#7AA6DCFF',Simple_repeat = "#003C67FF","Low_complexity" = "#3B3B3BFF")

colorManual <- c(LINE = '#A73030FF',SINE='#EFC000FF',LTR='#868686FF',"DNA"='#7AA6DCFF',Simple_repeat = "#003C67FF","Low_complexity" = "#3B3B3BFF")
otherAlpha <- 0.5
alphaManual <- c(LINE = 1,SINE=otherAlpha,LTR=otherAlpha,"DNA"=otherAlpha,Simple_repeat = otherAlpha,"Low_complexity" = otherAlpha)
x %>% 
	ggplot(aes(x = (`Number of peaks`), y = - `LogP enrichment (+values depleted)`,color = class,alpha = class)) +
	geom_point(size = 0.5) +
	geom_hline(yintercept = 2,color = "red",linetype = 2,alpha = 0.5, size =0.5)+
	ggrepel::geom_text_repel(data = label_combine,aes(label = Annotation ),size = 6/ggplot2::.pt,show.legend = F,segment.size = 0.1) +
	facet_wrap(~sample,scales = "free")+
	#geom_hline(yintercept = 0,color = "red",linetype = 2,alpha = 0.5, size =1)+
	#geom_vline(xintercept = 0,color = "red",alpha = 0.7,linetype = 2)+
	scale_color_manual(values = colorManual)+
	scale_alpha_manual(values = alphaManual)+
	#scale_size(range = c(0.1,2))+
	scale_y_log10()+
	xlab("Number of peaks")+ylab("Enrichment (-LogP value)")+
	guides(colour = guide_legend(override.aes = list(alpha = 1)))+
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
	)
ggsave("xaxis_numberPeaks_enrichment_high.pdf",width = 2.0,height = 3,useDingbats = F)

