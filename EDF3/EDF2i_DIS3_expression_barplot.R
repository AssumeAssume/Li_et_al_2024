library(tidyverse)
library(magrittr)
library(data.table)
library(glue)
library(DESeq2)
library(furrr)
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

setwd("/analysis2/lixf/RNA/human/HEK293/DIS3_muttion/Szczepinska_GenomeRes_2015_PRJNA270769/")

file <- dir("/analysis2/lixf/RNA/human/HEK293/DIS3_muttion/Szczepinska_GenomeRes_2015_PRJNA270769/FLL1/","t.test",full.names = T)

# REMOVE2 <- c(					"XRCC6","RBM15",	"ELOF1","INTS10","LEO1","RBPJ","ACTR8",
# 									"EWSR1","LSM10","HNRNPF","C7orf26","YEATS4","DDX42","NUP188","ZNHIT1","SRRD","SCAP")
df <- map(file,~fread(.x)) %>% reduce(rbind)
combineDf <- df 
combineDf <- combineDf %>% add_row(sample = "ctrl",FC=1,pvalue = 1)
combineDf %>% 
	mutate(sample = factor(sample, c("ctrl","RNB","PIN","PIN_RNB"))) %>% 
	ggplot()+
	geom_col(aes(x = sample, y = FC),width = 0.5,fill = "#743481")+
	geom_hline(yintercept = 1,color = "red",alpha = 0.7,linetype = 2,size = 0.3)+
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
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt")
	)+xlab("")+ylab("FoldChange")
ggsave("barplot_foldchange.pdf",width = 2,height  =1)
