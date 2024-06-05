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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/NG_revision/multiple_mapped_reads_individual")
LINE1_factor <- fread("/LiuLab/reference/Human/GRCh38/TE/LINE_1_factor.txt",header = F) %>% unlist
unique <- fread("unique_count.rawcounts.tsv")
multiple <- fread("multiple_count.rawcounts.tsv")
df <- left_join(unique,multiple,by='Geneid') 
colnames(df) <- c("uniqueID","unique","all")
plotdf <- df	%>% mutate(multiple = all-unique) %>% 
	filter(all > 5) %>% 
	mutate(subfamily = str_remove(uniqueID,"_chr.*")) %>% 
	mutate(UniquePercent = unique/all) %>% 
	mutate(subfamily = factor(subfamily,LINE1_factor) ) 
# plotdf %>% group_by(subfamily) %>% count() %>% fwrite("number_of_elements.tsv",sep = '\t')
all <- fread("./unique_count.rawcounts.tsv") %>% 
	mutate(subfamily =str_remove(Geneid,"_chr.*")) %>% group_by(subfamily) %>% 
	summarize(sum = n())
all
countdf <- plotdf %>% group_by(subfamily) %>% count() 
countdf %>% left_join(all) %>% 
	mutate(percent = n/sum) %>% 
	mutate(NumberOfElementsGT5 = n,
				 TotalNumber = sum) %>% 
		ggplot(aes(x = TotalNumber,y = NumberOfElementsGT5,color = percent,size = NumberOfElementsGT5))+
		geom_point()+
	scale_x_log10()+
	scale_y_log10()+
	scale_size(range = c(.1,2))+
	scale_colour_gradient2(mid = muted("blue"),high = muted("red"))+
	ggrepel::geom_text_repel(aes(label = subfamily),size = 6/ggplot2::.pt,
													 segment.size = 0.1,min.segment.length = 0)+
		theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 6,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "right",
				legend.key.width= unit(8, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt"),
				axis.line =  element_line(colour = "black", size = 0.25)
	)

ggsave("Number_of_elements_scatter.pdf",width = 3.5,height = 2.2,useDingbats= F)

