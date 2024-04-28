library(tidyverse)
library(magrittr)
library(data.table)
library(glue)
library(rstatix)
theme_set(theme_bw())
rm(list=ls())
filter <- dplyr::filter
select <- dplyr::select
desc <- dplyr::desc
reduce <- purrr::reduce
rename <- dplyr::rename
count <- dplyr::count
summarize <- dplyr::summarise
readr::local_edition(1)
options(stringsAsFactors = FALSE)
set.seed(629)
source("/analysis/lixf/my_script/my_R_function.R")

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure1_enhancer_transfection")
file <- dir("./","293T.txt")
names(file) <- c("293T")
df <- imap(file,function(x,y){
	tmp <- fread(x) %>% mutate(cellLine = y) 
	
	colnames(tmp) <- c("Control"   ,  "5UTR"      ,  "Down_5UTR"  , "5UTR(R)"   ,   "Down_5UTR(R)" ,"Oct4Enh"  ,     "Down-Oct4Enh",  "cellLine" )
	tmp <- tmp%>% 
		melt()
	tmp 
}) %>% reduce(rbind)
plotdf <- df %>% 
	mutate(cellLine = factor(cellLine,c("293T"))) %>%
	mutate(variable = factor(rev(variable),levels = rev(c("Control","5UTR","5UTR(R)","Oct4Enh","Down_5UTR","Down_5UTR(R)","Down-Oct4Enh")) )
	) %>%
	mutate(variable = rev(variable)) 
summarizeDf <-  plotdf %>% group_by(cellLine,variable) %>% 
	dplyr::summarize(    N             = length(value),
											 sd            = sd     (value),
											 se            = sd / sqrt(N),
											 mean          = mean(value),
											 ci            = se * qt(.95/2 + .5, N - 1)) %>% 	ungroup

ggplot()+
	geom_col(data = summarizeDf,aes(x = mean, y = variable),fill = "#808080",width = 0.6)+
	geom_jitter(data = plotdf,aes(x = value, y = variable),size = 0.9,stroke = 0.1,width = 0.1,height = 0.1)+
	geom_errorbar(data =summarizeDf, aes( x = mean, y = variable,xmin = mean - se,xmax = mean + se),colour="black", width=.4)+
	facet_wrap(~cellLine,scales = "free")+
	theme_classic()+
	scale_x_continuous(expand= c(0,0))+
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
	)+xlab("Relative activity") +ylab("")
ggsave("plot_enhancer_transfection293T.pdf",width = 4,height  = 2.7,useDingbats = F)

plotdf %>% filter(variable == "Control") %>% 
	select(cellLine,value)
plotdf %>% 
	group_by(cellLine) %>% 
	t_test()
plotdf %>% melt()
plotdf %>% 
	t_test(value ~ variable,p.adjust.method = "none") %>% 
	filter(str_detect(group2,"Control")) 
