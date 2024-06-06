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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure2/DIS3//")
l1 <- fread("/analysis/lixf/RNA/human/K562/KnockDown/DIS3/6-LINE1M/Combine_RepeatMasker_Unique.Deseq2.DIS3_vs_SAFE.p_0.05.host.annotated.bed")

l1_loc <- l1 %>% filter(str_detect(TE_ID,"^L1|^HAL1"))
l1_loc$Function %>% table()
tmp <- l1_loc %>% 
	mutate(localization = case_when(str_detect(Function,"ncRNA_intronic") ~ "intronic",
																	str_detect(Function,"ncRNA") ~ "ncRNA",
																	str_detect(Function,"stream") ~ "intergenic",
																	str_detect(Function,"UTR|splicing|exonic") ~ "protein coding gene",
																	TRUE ~ Function
	))
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
data$category
value = c("intronic" = "#435484","intergenic" = "#C95D3F","protein coding gene" = "#599583","ncRNA" = "#7BB8D2")
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
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
ggsave("DIS3_localization.pdf",width = 2,height = 1)
