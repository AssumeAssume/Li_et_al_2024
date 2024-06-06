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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/STARR/")
LINE1_factor <- fread("/LiuLab/reference/Human/GRCh38/TE/LINE_1_factor.txt",header = F) %>% unlist

file <- dir("./Repeatmasker_subfamily/","5UTR.tsv")
file <- file[!str_detect(file,"ENCFF")]
names(file) <- str_remove_all(file,"_random.*|.sorted.*")
df <- imap(file,function(x,y){
	tmp <- fread(glue("./Repeatmasker_subfamily/{x}"))
	colnames(tmp)[2] <- y
	return(tmp)
}) %>% reduce(left_join)

file <- dir("./RS_cCRE_subfamily//","5UTR.tsv")
file <- file[!str_detect(file,"ENCFF")]
names(file) <- str_remove_all(file,"_random.*|.sorted.*")
df2 <- imap(file,function(x,y){
	tmp <- fread(glue("./RS_cCRE_subfamily/{x}"))
	colnames(tmp)[2] <- y
	return(tmp)
}) %>% reduce(left_join)
df2_filter <- df2 %>% filter(!str_detect(Geneid,",CTCF-bound")) %>% 
	filter(str_detect(Geneid,"dELS|random")) #%>% 
#mutate(Geneid = ifelse(str_detect(Geneid,"ELS"),"ELS",Geneid)) %>% 
#group_by(Geneid) %>% summarise_all(sum) %>% ungroup

dfmelt <- rbind(df,df2_filter) %>% melt() %>% rename(Sample = variable)
TotalMapped <- fread("/analysis/lixf/tracks/STARR/several_cellLine/rawdata_analysis/0-QC/multiqc_data_3/multiqc_general_stats.txt") %>% 
	select(Sample,`Samtools_mqc-generalstats-samtools-mapped_passed`) %>% 
	rename(TotalMappedReads = `Samtools_mqc-generalstats-samtools-mapped_passed`) %>% 
	mutate(Sample = str_remove(Sample,"_random"))
plotdf <- dfmelt %>% left_join(TotalMapped) %>% 
	mutate(CPM = value*1e6/TotalMappedReads) %>% 
	mutate(Sample = case_when(str_detect(Sample,"K562_STARR_rep1_1") ~ "STARR_rep1 (PE100)",
														str_detect(Sample,"K562_STARR_rep1_2") ~ "STARR_rep1 (PE150)",
														str_detect(Sample,"K562_STARR_rep2") ~ "STARR_rep2 (PE100)",
														str_detect(Sample,"Ctrl_K562_HepG2") ~ "Control (PE150)",
														str_detect(Sample,"Combine_K562_STARR") ~ "STARR pooled",
	)) %>% 
	mutate(Sample = factor(Sample, levels = c("Control (PE150)", "STARR_rep1 (PE100)","STARR_rep2 (PE100)","STARR_rep1 (PE150)","STARR pooled"))) %>% 
	mutate(subfamily = factor(Geneid) ) %>% 
	mutate(subfamily = as.character(subfamily),
				 subfamily2 = ifelse(str_detect(subfamily,"L1PA8|L1PA[0-9]{2}"), "L1PA8-17",subfamily),
				 subfamily2 = factor(subfamily2, levels = c(LINE1_factor,"L1PA8-17","PLS","dELS","pELS","DNase-H3K4me3","DNase-only","Low-DNase","randomSequence"))
	)
p1 <- plotdf %>% 
	group_by(Sample,subfamily2) %>% 
	summarize(CPM = sum(CPM)) %>% 
	ggplot(aes(x = subfamily2 , y = CPM,color = Sample))+
	ggsci::scale_color_npg()+
	ylab("STARR-seq read counts (CPM)")+xlab("")+
	ggtitle("subfamily level")+
	geom_point(size = 0.5)+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "right",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt"),
				axis.line =  element_line(colour = "black", size = 0.25)
	)
p1
#ggsave(plot= p1,"./subfamily_STARR_CPM_l1pa8_17merge_withcCRE.pdf",width = 4.5,height = 2.4,useDingbats = F)

merge8_17 <- plotdf %>% 
	group_by(Sample,subfamily2) %>% 
	summarize(CPM = sum(CPM)) %>% ungroup
controlCPM <- merge8_17 %>% filter(Sample == "Control (PE150)") %>% select(subfamily2,CPM) %>% 
	rename(ControlCPM  =  CPM)

p3 <- merge8_17 %>% 
	#filter(Sample != "Control (PE150)") %>% 
	inner_join(controlCPM) %>% 
	mutate(FC = (CPM/ControlCPM)) %>% 
	ggplot(aes(x = subfamily2 , y = FC,color = Sample))+
	ggsci::scale_color_npg()+
	#geom_hline(yintercept = 0,color = "red",alpha = 0.7,linetype = 2)+
	ylab("Fold Change\n(compared with control)")+xlab("")+
	ggtitle("subfamily level")+
	geom_point(size = 0.5)+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "right",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt"),
				axis.line =  element_line(colour = "black", size = 0.25)
	)
p3

# individual --------------------------------------------------------------------------------------------------------------------------

file <- dir("./Repeatmasker_individual/","5UTR.tsv")
file <- file[!str_detect(file,"ENCFF")]
names(file) <- str_remove_all(file,"_random.*|.sorted.*")
df <- imap(file,function(x,y){
	tmp <- fread(glue("./Repeatmasker_individual/{x}"))
	colnames(tmp)[2] <- y
	return(tmp)
}) %>% reduce(left_join)

file <- dir("./RS_cCRE_individual//","5UTR.tsv")
file <- file[!str_detect(file,"ENCFF")]
names(file) <- str_remove_all(file,"_random.*|.sorted.*")
df2 <- imap(file,function(x,y){
	tmp <- fread(glue("./RS_cCRE_individual/{x}"))
	colnames(tmp)[2] <- y
	return(tmp)
}) %>% reduce(left_join)
df2
df2_filter <- df2 %>% filter(!str_detect(Geneid,",CTCF-bound")) %>% 
	filter(str_detect(Geneid,"dELS|random")) 
dfmelt <- rbind(df,df2_filter) %>% melt() %>% rename(Sample = variable)

TotalMapped <- fread("/analysis/lixf/tracks/STARR/several_cellLine/rawdata_analysis/0-QC/multiqc_data_3/multiqc_general_stats.txt") %>%
	select(Sample,`Samtools_mqc-generalstats-samtools-mapped_passed`) %>% 
	rename(TotalMappedReads = `Samtools_mqc-generalstats-samtools-mapped_passed`) %>% 
	mutate(Sample = str_remove(Sample,"_random"))
#plotdf$subfamily %>% table()


plotdf <- dfmelt %>% left_join(TotalMapped) %>% 
	filter(value > 5) %>% 
	mutate(CPM = value*1e6/TotalMappedReads) %>% 
	mutate(subfamily = str_remove_all(Geneid,"_chr.*|_EH.*")) %>% 
	#mutate(subfamily = ifelse(str_detect(subfamily,"ELS"),"ELS",subfamily)) %>% 
	mutate(Sample = case_when(str_detect(Sample,"K562_STARR_rep1_1") ~ "STARR_rep1 (PE100)",
														str_detect(Sample,"K562_STARR_rep1_2") ~ "STARR_rep1 (PE150)",
														str_detect(Sample,"K562_STARR_rep2") ~ "STARR_rep2 (PE100)",
														str_detect(Sample,"Ctrl_K562_HepG2") ~ "Control (PE150)",
														str_detect(Sample,"Combine_K562_STARR") ~ "STARR pooled",
	)) %>% 
	mutate(Sample = factor(Sample, levels = c("Control (PE150)", "STARR_rep1 (PE100)","STARR_rep2 (PE100)","STARR_rep1 (PE150)","STARR pooled"))) %>% 
	mutate(subfamily = factor(subfamily) ) %>% 
	mutate(subfamily = as.character(subfamily),
				 subfamily2 = ifelse(str_detect(subfamily,"L1PA8|L1PA[0-9]{2}"), "L1PA8-17",subfamily),
				 subfamily2 = factor(subfamily2, levels = c(LINE1_factor,"L1PA8-17","PLS","dELS","pELS","DNase-H3K4me3","DNase-only","Low-DNase","randomSequence"))
	)
#plotdf %>% select(subfamily,Geneid) %>% distinct() %>% group_by(subfamily) %>% count()
#plotdf %>% select(Geneid,subfamily2) %>% distinct() %>% group_by(subfamily2) %>% count()
p2 <- plotdf %>% 
	ggplot(aes(x = subfamily2 , y = CPM,fill = Sample))+
	ggsci::scale_fill_npg()+
	ylab("STARR-seq read counts (CPM)")+xlab("")+
	ggtitle("individual level")+
	geom_boxplot(width = 0.8,outlier.alpha = 0,lwd = 0.4)+
	coord_cartesian(ylim = c(0,1.0))+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "right",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt"),
				axis.line =  element_line(colour = "black", size = 0.25)
	)
p2

controlCPM <- plotdf %>% filter(Sample == "Control (PE150)") %>% select(Geneid,CPM) %>% 
	rename(ControlCPM  =  CPM)
p4 <- plotdf %>% 
	#filter(Sample != "Control (PE150)") %>% 
	inner_join(controlCPM) %>% 
	mutate(FC = (CPM/ControlCPM)) %>% 
	ggplot(aes(x = subfamily2 , y = FC,fill = Sample))+
	#geom_hline(yintercept = 1,color = "red",alpha = 0.7,linetype = 2)+
	ggsci::scale_fill_npg()+
	ylab("Fold Change\n(compared with control)")+xlab("")+
	ggtitle("individual level")+
	geom_boxplot(width = 0.8,outlier.alpha = 0,lwd = 0.4)+
	coord_cartesian(ylim = c(0,3))+
	#facet_wrap(~Sample)+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "right",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt"),
				axis.line =  element_line(colour = "black", size = 0.25)
	)
p4

((p1/p2)|(p3/p4))+plot_layout(guides = "collect")
ggsave("Combine_only_dELS_RS.pdf",width = 10,height = 5,useDingbats =F)

