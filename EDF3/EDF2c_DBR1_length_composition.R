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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure2/DBR1//")

background <- fread("/LiuLab/reference/Human/GRCh38/TE/Repeat_Masker_hg38_sorted.bed",header = F)
colnames(background) <- c("chr","start","end","strand","subfamily","class","family") 
back <- background %>% mutate(length = end-start) %>% mutate(group = rep("background",length(length)))
back_length_distribution <- back %>% select(subfamily:group)
back_L1 <- back_length_distribution %>% filter(family == "L1") %>% select(length,group)
de_l1 <- fread("/analysis/lixf/RNA/human/K562/KnockDown/DBR1/6-LINE1M/Combine_RepeatMasker_Unique.Deseq2.DBR1_vs_SAFE.p_0.05.L1.annotated.bed", fill=TRUE,header = T)
de_l1_df <- de_l1 %>% mutate(length = end - start,group = "sgDBR1 DE-L1", subfamily  = str_split_fixed(TE_ID,"_chr",2)[,1])

lengthdf <- rbind(back_L1,de_l1_df %>% select(length,group))
lengthdf %>% 
  ggplot() +
  geom_density(aes(x = length,group = group, color = group),size = 0.3)+
  ggsci::scale_color_aaas()+theme_bw()+
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
  )

ggsave("DBR1_length_distribution.pdf",height=1.5,width = 3,useDingbats = F)

