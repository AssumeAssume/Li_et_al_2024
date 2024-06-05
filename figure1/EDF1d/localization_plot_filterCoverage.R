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
readr::local_edition(1)
options(stringsAsFactors = FALSE)
set.seed(629)
source("/analysis/lixf/my_script/my_R_function.R")

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/K562_WGS_INSERTION/4-annoPlot")
anno <- fread("./resolution_1bp.L1insertion.removeContig_chrY.anno.txt") %>% mutate(Annotation  = str_remove_all(Annotation," \\(.*")) %>% 
	mutate(PeakID=paste(Chr,Start,End,sep ='_'))
bed <- fread("./resolution_1bp.L1insertion.removeContig_chrY.bed") %>% mutate(PeakID=paste(V1,V2+1,V3,sep ='_')) %>% 
	rename(Coverage = V4) %>% 
	select(PeakID,Coverage)

CoverageAnno <- anno %>% left_join(bed) %>% filter()
CoverageAnno %>% filter(Coverage>3) %>% fwrite("anno_filter_cov3.tsv",sep = '\t')


p <- plot_localization(CoverageAnno,3)

ggsave(plot = p, "filter_localization.pdf",width = 6,height = 3)

plot_localization <- function(df,n){
tmp <- df %>% filter(Coverage > n) %>% 
	mutate(localization = case_when(str_detect(Annotation,"exon|TTS") & `Gene Type` == "protein-coding" ~ "protein coding gene",
																	str_detect(Annotation,"exon|TTS") & `Gene Type` != "protein-coding" ~"ncRNA",
																	str_detect(Annotation,"intron")~"intronic",
																	str_detect(Annotation,"Intergenic")~"intergenic",
	)
	)

tmp2 <- tmp$localization %>% table() %>% as.data.frame()
colnames(tmp2) <- c("category","count")
data <- tmp2 %>% 
	mutate(category = fct_reorder(category,count,mean)) %>% 
	arrange(count)

# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$category, "\n N = ", data$count)
library(ggsci)
# Make the plot
value = c("intronic" = "#435484","intergenic" = "#C95D3F","protein coding gene" = "#599583","ncRNA" = "#7BB8D2")
data$category
p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
	geom_rect() +
	geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=6/ggplot2::.pt) + # x here controls label position (inner / outer)
	scale_fill_manual(values = value)+
	scale_color_manual(values = value) +
	coord_polar(theta="y") +
	xlim(c(-1, 4)) +
	theme_void() +
	theme(
		text=element_text(size= 8,family = "ArialMT",color = "black"),
	)+
	theme(legend.position = "none")
ggsave(plot = p,glue("localization_donuts_filter_{n}_Coverage.pdf"),width = 2,height = 1)
return(p)
}
