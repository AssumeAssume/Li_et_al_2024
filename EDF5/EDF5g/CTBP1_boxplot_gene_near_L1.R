#!/analysis2/software/miniconda3/bin/Rscript
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
readr::local_edition(1)
options(stringsAsFactors = FALSE)
set.seed(629)
source("/analysis/lixf/my_script/my_R_function.R")

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/Figure2h_RNA_K27ac")
DESeq2.knowngene.CTBP1_cKO <- fread("/analysis/lixf/RNA/human/K562/KO/CTBP1/5-DESeq2/KnownGene/Combine_KnownGene.Deseq2.CTBP1_vs_SAFE.result.txt") 

##### 
fusionGeneName <- fread("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure4_CTBP1/Promoter_L1_fusion.K27ac_increased_gene.name.txt",header = F) 
UP_L1 <- fread("./7-windowBed/up_p.0.05.L1/CTBP1_expression.tsv") %>% mutate(group = "RNA_up_L1")
file <- system("tree -f -i ./shuf_Inactive_FLL1P/7-windowBed/|grep -i CTBP1_expression.tsv",intern = T)
file
names(file) <- str_split_fixed(file,"/",10)[,4]
inactive <- imap_dfr(file,~fread(.x) %>% mutate(group = .y))


FC_df <- rbind(UP_L1,inactive)

distanceCutoff <- 32
plotdf <- FC_df %>% drop_na() %>%
	filter(! rowname %in% fusionGeneName$V1)  %>%
	mutate(distance = ifelse(abs(distance)<=distanceCutoff,0,distance)) %>% 
	mutate(distance = factor(distance)) %>% distinct()

plotdf %>% filter(group == "RNA_up_L1")%>% group_by(distance) %>% dplyr::count() %>% fwrite("CTBP1KO.count.number.tsv",sep ='\t' )

MedianInactive <- plotdf %>% 
	filter(str_detect(group,"shuf")) %>% 
	mutate(group2 = str_remove(group,"_[0-9]+")) %>% 
	group_by(group2,distance,group) %>% 
	summarize(Q3 = quantile(log2FoldChange,0.75),
						Q1 = quantile(log2FoldChange,0.25),
						Median = median(log2FoldChange)
						) %>% 
	ungroup %>% 
	group_by(distance,group2) %>% 
	summarize(MedianQ3 = median(Q3),
						MedianQ1 = median(Q1),
						MedianMedian = median(Median)
	) %>% ungroup

plotdf %>% 
	filter(str_detect(group,"RNA_up_L1")) %>% 
	ggplot(aes(x = distance, y = log2FoldChange))+
	geom_boxplot(outlier.alpha = 0,lwd =0.3)+
	ungeviz::geom_hpline(data= MedianInactive,aes(y = MedianMedian),colour = "red",width = 0.5,size=0.6,)+
	ungeviz::geom_hpline(data= MedianInactive,aes(y = MedianQ3),colour = "red",width = 0.2,size=0.6,alpha=0.5)+
	ungeviz::geom_hpline(data= MedianInactive,aes(y = MedianQ1),colour = "red",width = 0.2,size=0.6,alpha = 0.5)+
	geom_segment(data= MedianInactive,aes(x=distance,y = MedianMedian,xend = distance,yend = MedianQ3),color="red",size = 0.3,alpha = 0.5)+
	geom_segment(data= MedianInactive,aes(x=distance,y = MedianMedian,xend = distance,yend = MedianQ1),color="red",size = 0.3,alpha = 0.5)+
	#geom_hline(yintercept = 0,color = "black",alpha = 0.7,linetype = 1)+
	coord_cartesian(ylim = c(-1.55,1.9))+
	scale_y_continuous(expand = c(0,0))+
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
	)+guides(fill=guide_legend(title=""))  + xlab("Distance to L1 [kb]")
ggsave("RNA_up_L1_with_Inactive_FLL1P.pdf",width = 3,height = 2)
