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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/K562_WGS_INSERTION/3-Sam2Bam")
tmp <- fread("tmp2",header = F)
tmp3 <- tmp %>% 
	mutate(left =  case_when(str_detect(V7,"Left") ~ str_replace_all(V7,"Left([0-9]+)S.*","\\1"),
													 TRUE ~ "0"
													 ) %>% as.numeric,
				 Right =  case_when(str_detect(V7,"Right") ~ str_replace_all(V7,".*Right([0-9]+)S","\\1"),
				 									TRUE ~ "0") %>% as.numeric
				 ) %>% 
	mutate(delta = left-Right)

tmp3 %>% 
	group_by(V12) %>% 
	mutate(sum = sum(delta)) %>% 
	mutate(V9 = ifelse(sum < 0,V10-1,V9),
				 V10 = ifelse(sum > 0,V9+1,V10)
				 ) %>% 
	select(V8:V12) %>% 
	distinct() %>% 
	fwrite("resolution_1bp.L1insertion.bed",sep ='\t',col.names = F)
