library(tidyverse)
library(magrittr)
library(data.table)
library(glue)
library(ungeviz) # geom_hpline
library(patchwork)
theme_set(theme_bw())
rm(list=ls())
filter <- dplyr::filter
select <- dplyr::select
desc <- dplyr::desc
reduce <- purrr::reduce
rename <- dplyr::rename
options(stringsAsFactors = FALSE)
set.seed(629)
source("/analysis/lixf/my_script/my_R_function.R")

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure1/figure1_20220805/")

act_fdr_FACS <- fread("../LINE1_activator_fdr_0.1_FACS_20220803.csv")%>% mutate(group = "supporter") 
colnames(act_fdr_FACS)[1] <- "id"
act_fdr_FACS <- mutate(act_fdr_FACS,id2=fct_reorder(id,desc(value),mean)) 
act_fdr_FACS$id %>% table() %>% table

sup_fdr_FACS <- fread("../LINE1_suppressor_fdr_0.1_FACS_20220803.csv") %>% mutate(group = "suppressor") 
colnames(sup_fdr_FACS)[1] <- "id"
sup_fdr_FACS <- mutate(sup_fdr_FACS,id2=fct_reorder(id,desc(value),mean)) 

sup_fdr_FACS$id %>% table() %>% table()

df <- rbind(act_fdr_FACS,sup_fdr_FACS)
df$id2 <- factor(df$id2,c(levels(sup_fdr_FACS$id2), levels(act_fdr_FACS$id2)) %>% unique)

filterDf <- df
#### add another group of NC
library(rstatix)
tmp <- filterDf %>% 
	mutate(group = factor(group,levels = c("suppressor","supporter"))) %>% 
	filter(group == "suppressor") %>% 
	filter(id != "NC") %>% 
	select(id2,value) %>% mutate(group = "KO")
length <- tmp$id2 %>% unique %>% length()
ncValue <- tmp %>% filter(id == "NC") %>% pull(value) %>% unique()
NC <- data.frame(id2 = rep(tmp$id2 %>% unique, each =  3 ),value = rep(c(ncValue),length)) %>% mutate(group = "NC")
suppressor_pvalue <- rbind(tmp,NC) %>% 
	group_by(id2) %>% 
	t_test(  value ~ group,alternative = "greater")

fwrite(suppressor_pvalue,"suppressor_FACS_pvalue_df.tsv",sep = '\t')

colorManual <- c(supporter = "#b92b27",suppressor = "#1565C0")
Psupp <- filterDf %>% 
	mutate(group = factor(group,levels = c("suppressor","supporter"))) %>% 
	filter(id != "NC") %>% 
	filter(group == "suppressor") %>% 
	mutate(id2=fct_reorder(id,desc(value),mean)) %>%
	ggplot(aes(x=id2,y=value,color = group))+
	#geom_bar(stat="identity")+
	stat_summary(fun = "mean", colour = "black",width = 0.5,size=0.6,geom = "hpline",alpha=0.8)+
	geom_hline(yintercept = 1,linetype=2,color="black",alpha=0.5,size=0.2)+
	geom_point(shape=1,size=0.5)+
	scale_color_manual(values = colorManual)+
	#facet_wrap(~group,scales = "free",nrow = 2)+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.ticks.length.x = unit(0.05, "cm"),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "bottom",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
	)+guides(fill=guide_legend(title=""),color = "none")+xlab("")+ylab("FACS GFP value \n(relative to NC)")


library(rstatix)
tmp <- filterDf %>% 
	mutate(group = factor(group,levels = c("suppressor","supporter"))) %>% 
	filter(group == "supporter") %>% 
	filter(id != "NC") %>% 
	select(id2,value) %>% mutate(group = "KO")
length <- tmp$id2 %>% unique %>% length()
NC <- data.frame(id2 = rep(tmp$id2 %>% unique, each =  3 ),value = rep(c(ncValue),length)) %>% mutate(group = "NC")
supporter_pvalue <- rbind(tmp,NC) %>% 
	group_by(id2) %>% 
	t_test(  value ~ group,alternative = "less")

fwrite(supporter_pvalue,"activator_FACS_pvalue_df.tsv",sep = '\t')

Pacti <- filterDf %>% 
	mutate(group = factor(group,levels = c("suppressor","supporter"))) %>% 
	filter(id != "NC") %>% 
	filter(group == "supporter") %>% 
	mutate(id2=fct_reorder(id,desc(value),mean)) %>%
	ggplot(aes(x=id2,y=value,color = group))+
	#geom_bar(stat="identity")+
	stat_summary(fun = "mean", colour = "black",width = 0.5,size=0.6,geom = "hpline",alpha=0.8)+
	geom_hline(yintercept = 1,linetype=2,color="black",alpha=0.5,size=0.2)+
	geom_point(shape=1,size=0.5)+
	scale_color_manual(values = colorManual)+
	#facet_wrap(~group,scales = "free",nrow = 2)+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.ticks.length.x = unit(0.05, "cm"),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "bottom",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
	)+guides(fill=guide_legend(title=""),color = "none")+xlab("")+ylab("FACS GFP value \n(relative to NC)")

Psupp/Pacti
ggsave("GFP_FACS_combine_usingMean_230102.pdf",height=2.8,width = 3.2,useDingbats = F)
