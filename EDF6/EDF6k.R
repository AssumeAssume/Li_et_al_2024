library(tidyverse)
library(magrittr)
library(data.table)
library(glue)
library(scales)
library(patchwork)
theme_set(theme_bw())
rm(list=ls())
filter <- dplyr::filter
select <- dplyr::select
desc <- dplyr::desc
reduce <- purrr::reduce
options(stringsAsFactors = FALSE)
set.seed(629)
source("/analysis/lixf/my_script/my_R_function.R")

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure2/MTC")

############### exp/obs enrichment file ################### 
file <- "/analysis2/lixf/m6A/human/ChIP/Xu_MolecularCell_2022_m6A_MCF7_HEK293T/3-MACS2/MCF7_merge.narrowPeak_enrichment.txt"
annlist <- map(file,~ fread(.x,skip = 6))


####### filter out copies < 10, filter Promoter intron ..., only remain repeats, plot log2 Ratio (obs/exp)
family <- fread("/LiuLab/reference/Human/GRCh38/TE/Repeat_Masker_hg38.class",header = F)
colnames(family) <- c("Annotation","class","family")

filterlist <- map(annlist,~ .x %>% 
										right_join(family) %>% 
										filter(`Number of peaks` > 5 & 
													 	! Annotation == "ALR/Alpha" & 
													 	!str_detect(family,"\\?") & 
													 	!str_detect(class,"Satellite")
										)
)


x <- filterlist[[1]]
y <- "MTC"

x <- x %>% mutate(class = ifelse(class %in% c("LINE","SINE","LTR","Retroposon","DNA","SVA","Simple_repeat","Low_complexity"),class, "others"))

# enrichment number  ------------------------------------------------------------------------------------------------------------------------


label <- x %>% top_bottom(5,`Log2 Ratio (obs/exp)`)
label2 <- x %>% top_bottom(5,`LogP enrichment (+values depleted)`)
label_family <- x %>% group_by(class) %>% top_bottom(2,`Log2 Ratio (obs/exp)`) %>% ungroup
label_L1 <- x %>% filter(str_detect(Annotation ,"L1PA3"))
label_combine <- reduce(list(label_family,label,label2,label_L1),rbind) %>% distinct()
colorManual <- c(LINE = '#A73030FF',SINE='#EFC000FF',LTR='#0073C2FF',Retroposon='#868686FF',"DNA"='#7AA6DCFF',Simple_repeat = "#003C67FF","Low_complexity" = "#3B3B3BFF",others = "black")
x %>% 
	ggplot(aes(x = `Number of peaks`, y = - `LogP enrichment (+values depleted)`,color = class)) +
	geom_point(size = 0.5) +
	geom_hline(yintercept = 2,color = "red",linetype = 2,alpha = 0.5, size =1)+
	ggrepel::geom_text_repel(data = label_combine,aes(label = Annotation ),size = 6/ggplot2::.pt,show.legend = F) +
	#geom_hline(yintercept = 0,color = "red",linetype = 2,alpha = 0.5, size =1)+
	#geom_vline(xintercept = 0,color = "red",alpha = 0.7,linetype = 2)+
	ggtitle("MTC enrichment over repeats")+
	scale_color_manual(values = colorManual)+
	#scale_size(range = c(0.1,2))+
	scale_y_log10()+
	xlab("Number of peaks")+ylab("Enrichment (-LogP value)")+
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
	)

ggsave("MTC_peaks_enrichment.pdf",height=3,width = 4.5,useDingbats = F)
