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



setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/Referee3/Figure2g_add_down_gene/absolute_distance")


file <- dir("./background_exp","distance",full.names = T)


backgroundTmp <- imap(file, function(x,y){
	tmp <- fread(x) %>% filter(V13 != -1) %>% select(V13)  %>% 
		mutate(bin = cut_width(V13, 10000,boundary = 0,label =F)) 
	a <- tmp %>% filter(bin == 1) %>% head(1) %>% pull(V13)
	tmp <- tmp %>% mutate(bin = bin+a%/%10000)%>% 
		pull(bin) %>% table() 
	tmpp <- data.frame(bin = names(tmp) %>% as.numeric, Freq = tmp %>% as.numeric) %>% 
		filter(bin <=100)
	
	colnames(tmpp) <- c("bin",paste0("group",y))
	tmpp
}
) %>% reduce(full_join)
backgroundTmp
backgroundTmp[is.na(backgroundTmp)] <- 0
tmpDf <- backgroundTmp %>% arrange(bin) %>% filter(bin <=100) 
mat <- tmpDf[,-1]
rownames(mat) <- tmpDf$bin

conf.int <- apply(mat,1,function(x) t.test(x,conf.level = 0.95)$conf.int)
ConfBackground <- conf.int %>% t() %>% as.data.frame() %>% rownames_to_column(var =  "Var1")
fwrite(ConfBackground,"./down_background_Conf.tsv",sep ='\t' )



# up background conf ------------------------------------------------------------------------------------------------------------------



file <- dir("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure4_CTBP1/absolute_distance//background_exp","distance",full.names = T)


backgroundTmp <- imap(file, function(x,y){
	tmp <- fread(x) %>% filter(V13 != -1) %>% select(V13)  %>% 
		mutate(bin = cut_width(V13, 10000,boundary = 0,label =F)) 
	a <- tmp %>% filter(bin == 1) %>% head(1) %>% pull(V13)
	tmp <- tmp %>% mutate(bin = bin+a%/%10000)%>% 
		pull(bin) %>% table() 
	tmpp <- data.frame(bin = names(tmp) %>% as.numeric, Freq = tmp %>% as.numeric) %>% 
		filter(bin <=100)
	
	colnames(tmpp) <- c("bin",paste0("group",y))
	tmpp
}
) %>% reduce(full_join)
backgroundTmp
backgroundTmp[is.na(backgroundTmp)] <- 0
tmpDf <- backgroundTmp %>% arrange(bin) %>% filter(bin <=100) 
mat <- tmpDf[,-1]
rownames(mat) <- tmpDf$bin

conf.int <- apply(mat,1,function(x) t.test(x,conf.level = 0.95)$conf.int)
ConfBackground <- conf.int %>% t() %>% as.data.frame() %>% rownames_to_column(var =  "Var1")
fwrite(ConfBackground,"./up_background_Conf.tsv",sep ='\t' )
