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
plotdf %>% group_by(subfamily) %>% count() %>% fwrite("number_of_elements.tsv",sep = '\t')

countdf <- plotdf %>% group_by(subfamily) %>% count() 
TotalN <- countdf$n %>% reduce(sum)
plotdf %>% left_join(countdf) %>% mutate(subfamily = paste0(subfamily,"\n","(n = ",n,")"))
library(MetBrewer)
MetBrewer::colorblind_palettes
levelsDf <- data.frame(subfamily = levels(plotdf$subfamily))	%>% left_join(countdf) %>% drop_na() %>% 
	mutate(subfamily = paste0(subfamily,"\n","(n = ",n,")"))
plotdf %>% 
	left_join(countdf) %>% 
	mutate(subfamily = paste0(subfamily,"\n","(n = ",n,")"),
				 subfamily = factor(subfamily,levels = levelsDf$subfamily)) %>% 
	ungroup %>% 
	#mutate(variable =fct_inorder(variable) %>% fct_rev) %>% 
	ggplot(aes(x = subfamily , y = UniquePercent))+
	geom_boxplot(outlier.alpha = 0,width = 0.7,lwd = 0.25)+
	#ggbeeswarm::geom_quasirandom(size=0.3,stroke = 0.01)+
	scale_y_continuous(label = percent,breaks = seq(0,1,0.2))+
	xlab("")+ylab("Unique mapped read fraction (%)")+
	#scale_fill_manual(values = c("#1A1919","#EE0000FF"),name="")+
	#scale_fill_aaas()
	#scale_fill_manual(values = met.brewer("VanGogh3",7)[c(1,3)],name="") +
	theme_classic()+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 6,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "bottom",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt"),
				axis.line =  element_line(colour = "black", size = 0.25)
	)+ggtitle(glue("Total n = {TotalN}"))+labs(caption = "only elements with read counts \ngreater than 5 were kept for analysis")

ggsave("reads_L1_percent_withNumber.pdf",width = 3,height = 2,useDingbats= F)

