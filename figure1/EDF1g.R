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
options(stringsAsFactors = FALSE)
set.seed(629)
source("/analysis/lixf/my_script/my_R_function.R")
setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure1/figure1//")
sup_fdr <- L1_sup %>% filter(GFP_plus.fdr < 0.1 ) %>% select(id,GFP_plus.lfc)
sup_sgRNA <- fread("/analysis/lixf/sgRNA/human/GFP/results/test/rra.sgrna_summary.txt")
fdr_sup_sgRNA <- sup_sgRNA %>% filter(Gene %in% sup_fdr$id)
sup_sgRNA_combine <- fdr_sup_sgRNA %>% group_by(Gene) %>% summarise(mean = mean(LFC),
																																		median = median(LFC),
																																		max = max(LFC),
																																		min = min(LFC)) %>% arrange(desc(mean))


act_fdr <- L1_act %>% filter(GFP_minus.fdr < 0.1 ) %>% select(id,GFP_minus.lfc) %>% mutate(GFP_minus.lfc = -1*GFP_minus.lfc)
act_sgRNA <- fread("/analysis/lixf/sgRNA/human/GFP_2/results/test/rra.sgrna_summary.txt")
fdr_act_sgRNA <- act_sgRNA %>% filter(Gene %in% act_fdr$id)
act_sgRNA_combine <- fdr_act_sgRNA %>% group_by(Gene) %>% summarise(mean = -mean(LFC),
																																		median = -median(LFC),
																																		max = -max(LFC),
																																		min = -min(LFC)) %>% arrange((median))

combine <- rbind(sup_sgRNA_combine,act_sgRNA_combine)
colnames(combine) <- c("Gene","mean_lfc","median_lfc","max_lfc","min_lfc")
n = 80
star_label <- c("MORC2","FAM208A","SETDB1","SAFB","MPHOSPH8","ATF7IP","PPHLN1","YY1","MYC")
top <- combine %>%
	mutate(Gene = ifelse(Gene %in% star_label,paste0(Gene," *"),Gene)) %>% 
	ungroup %>% 
	top_bottom(n,mean_lfc) %>% 
	arrange(desc(mean_lfc)) %>% 
	mutate(Gene = fct_inorder(Gene)) 

label_df <- top
label_df <- label_df %>% mutate(label_position = NA) %>% relocate(label_position,.before = mean_lfc)
#colnames(label_df)[2] <- c("label_position")
distance <- 1.5
for (i in 1:nrow(top)) {
	if (top[i,3] > 0 ) {
		if(i %% 2 == 1 ){
			label_df[i,2] <- top[i,5]- distance
		}else {
			label_df[i,2] <- top[i,4]+ distance
		}
	}else if( top[i,2] <= 0 ){
		if(i %% 2 == 1 ){
			label_df[i,2] <- top[i,5]+ distance
		}else {
			label_df[i,2] <- top[i,4]- distance
		}
	}
}

top %>%
	ggplot(aes(x = Gene, y = mean_lfc))+
	geom_point(aes(color  = ifelse(mean_lfc > 0, "supp","acti"),fill = ifelse(mean_lfc > 0, "supp","acti")),shape = 23,size = 0.5)+
	geom_hline(yintercept = 0,alpha = 0.6,size = 0.8)+
	geom_vline(xintercept = n+0.5,alpha = 0.5,size = 0.3,linetype = 2)+
	geom_segment(aes(x =Gene, y = mean_lfc, xend = Gene, yend = max_lfc ,color = ifelse(mean_lfc > 0, "supp","acti")),size = 0.25)+
	geom_segment(aes(x =Gene, y = mean_lfc, xend = Gene, yend = min_lfc ,color = ifelse(mean_lfc > 0, "supp","acti")),size = 0.25)+
	geom_text(data = label_df,aes(x = Gene, y = label_position,label = Gene,color = ifelse(mean_lfc > 0, "supp","acti")),size = 6/ggplot2::.pt,angle = 90)+
	scale_y_continuous(limits = c(-6,6),breaks = c(-6,seq(-4,4,2),6),labels = abs(c(-6,seq(-4,4,2),6)),expand = c(0,0.01))+
	scale_color_manual(values = c("#b92b27","#1565C0"))+
	scale_fill_manual(values = c("#b92b27","#1565C0"))+
	scale_x_discrete(guide = guide_axis(n.dodge=1))+
	theme_classic()+
	labs(x = "GFP screen with L1-UTR-reporter in K562 cells",
			 y = expression(log["2"]~"Fold enrichment [ GFP-sort / Ctrl ]"))+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 6,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_blank(),
				axis.line.x.bottom = element_blank(),
				axis.ticks.x=element_blank(),
				axis.ticks.length.x = unit(0.05, "cm"),
				axis.ticks.x.top = element_line(size = 0.1),
				axis.line.x.top = element_line(size = 0.1 ),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "bottom",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
	)+
	guides(color=guide_legend(title="",override.aes = list(size=1)),fill="none",size="none")


ggsave("./top80_screen_lfc_rankplot_using_mean_Arial_6pt.pdf",width = 8.5,height = 3)
