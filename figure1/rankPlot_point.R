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


setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure1/figure1_20220805/")

# activator ---------------------------------------------------------------------------------------------------------------------------

df <- read_tsv("/analysis/lixf/sgRNA/human/GFP_2/results/test/rra.gene_summary.txt",col_names = T)
########## suppressor
genename2 <- c("METTL3","METTL14","ZC3H13","ACTL6A","YY1","RBM15","GTF3C4","FAM208B")

genename <- c("")

supp <- df %>%
	arrange(`pos|rank`) %>% 
	mutate(rank=`pos|rank`,
				 pvalue=ifelse(`pos|p-value`==0,minpvalue/5,`pos|p-value`)
	)

labelaa <- supp %>%filter(id %in% genename2)
labelred <- supp %>% filter(id %in% genename)
labeldf <- full_join(labelaa,labelred)

top_100_score <- (supp %>% head(100) %>% tail(1))$`pos|score` 
logtop_100_score <- -log(top_100_score,10)

fdrcutoff <- 0.1
cutoffline <- supp %>% filter(`pos|fdr` < fdrcutoff) %>% tail(1) %$%`pos|score` 

p2 <- supp %>% 
	ggplot(aes(x=rank,y=-log10(`pos|score` )))+
	geom_point(size=0.1,stroke = 0.1,
						 alpha = 0.7,color = "grey")+
	geom_point(data = supp %>% filter(id %in% genename2),size =0.5,stroke = 0.1,
						 alpha = 1,color ="black")+
	geom_hline(yintercept = -log10(cutoffline),color="black",linetype=2,size = 0.1)+
	geom_text(data=data.frame(x=13000,y=-log10(cutoffline)), aes(x, y), label=paste0("FDR cutoff: ",fdrcutoff), size=6/ggplot2::.pt,vjust=-1)+
	ggrepel::geom_text_repel(data=labeldf,label=labeldf$id,size=6/ggplot2::.pt, force = 1,segment.size=0.05,segment.alpha=0.8,segment.color=ifelse(labeldf$id%in%genename,"red","black"),color=ifelse(labeldf$id%in%genename,"red","black"),direction="both",min.segment.length = 0)+
	theme_classic()+
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
	)+
	scale_x_continuous(breaks = c(0,5000,10000,15000,20000))+
	xlab("Gene rank")+ylab("-log10( RRA score)")#ylim(c(0,6))
p2
pp <- p2&
	theme(
		axis.ticks = element_line(colour = "black", size = 0.1),
		axis.line = element_line(colour = 'black', size = 0.1)
	)


ggsave(plot= pp,"activator_rankplot_fdr_0.1.pdf",height=1.5,width = 2,useDingbats = F)


