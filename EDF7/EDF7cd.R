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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure5_NCCIT_CRISPRai")
total_FL <- fread("/LiuLab/reference/Human/GRCh38/TE/TE_bed/LINE-1.FullLength.strand.bed") %>% 
	mutate(uniqueID=paste(V4,V1,V2,V3,V6,sep = "_"))

K27ac_up_L1 <- fread("/analysis/lixf/CRISPRa_i/NCCIT/ChIP/H3K27ac/CRISPRa/4-overlapRepeats/H3K27ac/sgL1_up/up_L1.bed") %>% 
	mutate(uniqueID=paste(V4,V1,V2,V3,V6,sep = "_"),
				 lengthGroup = ifelse(V3-V2 >=6000, "FL","non-FL")) %>% 
	filter(lengthGroup == "FL")

mismatch <- c(0,1,2,3)
m <- 0

list <- map(mismatch,function(m){
	
dir <- glue("/analysis/lixf/cas9_design/L1_dcas9_target_cal/4-OnTarget/mismatch_{m}")
file <- dir(path = dir,"HS-[578].*ontarget.candidates.bed$",full.names = T)
target_L1 <- map(file,~ fread(.x,header = F) %>% select(-1:-12)) %>%  
	reduce(rbind) %>% distinct() %>%
	mutate(uniqueID=paste(V16,V13,V14,V15,V18,sep = "_"),
				 lengthGroup = ifelse(V15-V14 >=6000, "FL","non-FL")) %>%
	filter(lengthGroup == "FL")

target_L1 %>% filter(str_detect(V16,"L1HS|L1PA[23]"))

#RNA_up_L1 <- fread("/analysis/lixf/CRISPRa_i/NCCIT/RNA/WY_20220719/CRISPRa/6-LINE1M/up_L1.bed")

k27ac_target_fll1 <- target_L1 %>% mutate(group = ifelse(uniqueID %in% K27ac_up_L1$uniqueID,"K27ac + target","K27ac - target")) %>%
	select(uniqueID,group)
combine_df <- left_join(total_FL,k27ac_target_fll1)
combine_df[is.na(combine_df)] <- "Nontarget"
LINE1_factor <- fread("/LiuLab/reference/Human/GRCh38//TE/UCSCrmsk/L1.div.txt",header = F) %>% pull(V1) %>% unlist()
# combine_df %>% 
# 	mutate()

p <- combine_df %>%
	group_by(V4,group) %>%
	dplyr::count() %>%
	ungroup %>%
	filter(str_detect(V4,"L1HS|L1PA[2-8]|L1MA[1-3]")) %>%
	mutate(V4 = factor(V4,levels = rev(LINE1_factor) ),
				 group = factor(group,levels = rev(c("K27ac + target","K27ac - target","Nontarget")) )
				 ) %>%
	ggplot(aes(x = V4, y = n,fill = group))+
	geom_bar(stat = "identity",position = "stack") +
	coord_flip()+
	scale_fill_manual(values = rev(c("#0073C2FF" ,"#EFC000FF", "#868686FF")),name = "")+
	ggtitle(glue("mismatch {m}"))+
	labs(y = "Counts",x= "")+
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
return(p)
})
mypal <- ggsci::pal_jco()(3)
show_col(mypal)

reduce(list,`|`)+plot_layout(guides = "collect")
ggsave("target_calculation_K27ac.pdf",height=2,width = 8,useDingbats = F)


# without K27ac -----------------------------------------------------------------------------------------------------------------------
mismatch <- c(0,1,2,3)
list <- map(mismatch,function(m){
	dir <- glue("/analysis/lixf/cas9_design/L1_dcas9_target_cal/4-OnTarget/mismatch_{m}")
	file <- dir(path = dir,"HS-[578].*ontarget.candidates.bed$",full.names = T)
	target_L1 <- map(file,~ fread(.x,header = F) %>% select(-1:-12)) %>%  
		reduce(rbind) %>% distinct() %>%
		mutate(uniqueID=paste(V16,V13,V14,V15,V18,sep = "_"),
					 lengthGroup = ifelse(V15-V14 >=6000, "FL","non-FL")) %>%
		filter(lengthGroup == "FL")
	
	
	#RNA_up_L1 <- fread("/analysis/lixf/CRISPRa_i/NCCIT/RNA/WY_20220719/CRISPRa/6-LINE1M/up_L1.bed")
	
	k27ac_target_fll1 <- target_L1 %>% mutate(group = "Target") %>%
		select(uniqueID,group)
	combine_df <- full_join(total_FL,k27ac_target_fll1)
	combine_df[is.na(combine_df)] <- "Nontarget"
	LINE1_factor <- fread("/LiuLab/reference/Human/GRCh38//TE/UCSCrmsk/L1.div.txt",header = F) %>% pull(V1) %>% unlist()
	p <- combine_df %>%
		group_by(V4,group) %>%
		count() %>%
		ungroup %>%
		filter(str_detect(V4,"L1HS|L1PA[2-8]|L1MA[1-3]")) %>%
		mutate(V4 = factor(V4,levels = rev(LINE1_factor) ),
					 group = factor(group,levels = rev(c("Target","Nontarget")) )
		) %>%
		ggplot(aes(x = V4, y = n,fill = group))+
		geom_bar(stat = "identity",position = "stack") +
		coord_flip()+
		scale_fill_manual(values = rev(c("#0073C2FF" , "#868686FF")),name = "")+
		ggtitle(glue("mismatch {m}"))+
		labs(y = "Counts",x= "")+
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
	return(p)
})

reduce(list,`|`)+plot_layout(guides = "collect")
ggsave("target_calculation_withoutK27ac.pdf",height=2,width = 8,useDingbats = F)


