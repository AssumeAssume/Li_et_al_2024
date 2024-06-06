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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure2/KMT2D/")
R <- fread("/analysis2/lixf/RNA/mouse/Liver_tumor/Kmt2d/Wang_CancerDiscovery_2020_Kmt2d_mouse/5-DESeq2/RepeatMasker/cancer_discovery.tsv")
R <- R %>% 
	mutate(class = str_split_fixed(gene,":",3)[,2],
				 rowname = str_split_fixed(gene,":",3)[,1])
labeldf <- R %>% filter(class =="L1" & padj < 0.001 )


LINE1_factor <- read.table("/LiuLab/reference/Mouse/GRCm38/TE/LINE_1_factor.txt",col.names = F) %>% as.list() %>% unlist
L1_select <- c("L1Md|L1_Mus|L1_Mur|L1VL")
plotdf <- R %>% 
	filter( str_detect(rowname,"^L1")) %>% 
	filter(str_detect(rowname,L1_select)) %>% 
	mutate(LINE1=factor(rowname,levels=LINE1_factor)) %>% 
	mutate(pvalue_group = case_when(pvalue <= 0.01 ~ "<= 0.01",
																	pvalue > 0.01 & pvalue <= 0.05  ~ "<= 0.05",
																	TRUE ~ "> 0.05"),
				 sample = "Kmt2d KO MAL1 cell"
	) %>% 
	arrange(desc(log2FoldChange)) 

sizeManual <- c("<= 0.01" = 3, "<= 0.05" = 2.3 , "> 0.05" = 1.0)
max <- plotdf$log2FoldChange %>% max()
min <- plotdf$log2FoldChange %>% min()
range <- max - min
midPoint <- 1 - max/range


#### LINE-1 ------
p <- plotdf %>%
	ggplot(aes(x=LINE1, y=sample, size=pvalue_group, fill=log2FoldChange)) +
	geom_point(alpha=1, shape=21, color="black",stroke=0.005) +
	scale_size_manual(values = sizeManual ,name="pvalue") +
	scale_fill_gradient2( low = "darkblue",mid = "white",high = muted("red"),
												breaks=c(-1,-0.2,0,0.2,0.4,0.6),limits = c(-0.25,0.8)
												)+
	# theme_minimal() +
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "bottom",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt")
	)+
	guides(fill = guide_colorbar(title.position = "top",
															 title.hjust = .5,
															 barwidth = unit(5, "lines"),
															 barheight = unit(.5, "lines")),
				 size = guide_legend(title.position = "top",
				 										title.hjust = .5,
				 										#barwidth = unit(10, "lines"),
				 										#barheight = unit(.5, "lines")
				 )
	)+
	xlab("")+ylab("") 
p
ggsave(plot = p,"./selected_pvalue_LINE1_foldchange.pdf",width = 3.5,height =1.6,useDingbats=FALSE)  
#guides(fill=guide_legend(title="padj"))+

