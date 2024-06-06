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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/figure4_CTBP1")

R <- fread("/analysis/lixf/RNA/human/K562/KO/CTBP1/5-DESeq2/RepeatMasker_Unique/Combine_RepeatMasker_Unique.Deseq2.CTBP1_vs_SAFE.result.txt")
gene <- fread("/analysis/lixf/RNA/human/K562/KO/CTBP1/5-DESeq2/KnownGene/Combine_KnownGene.Deseq2.CTBP1_vs_SAFE.result.txt")
annotation <- vroom::vroom("/LiuLab/reference/Human/GRCh38/TE/TE_annotation.tsv",col_names = T)

DESeqdf <- R

	df <- as.data.frame(DESeqdf)
	df$padj[is.na(df$padj)] <- 1
	colnames(df)[1] <- "symbol"
	data <- df 
	
	repeatclass <- left_join(data %>% mutate(anno=symbol),
													 annotation %>% dplyr::select(anno,X5,X6,X7),
													 by="anno")
	tmp <- repeatclass %>% mutate(repeatclass=X7,
																repeatclass = case_when(str_detect(anno,"L1HS|L1PA2|L1PA3") ~ "L1PA1-3",
																												str_detect(anno,"L1PA") ~ "other L1PA",
																												repeatclass == "L1" ~ "nonL1PA LINE1",
																												TRUE ~ repeatclass
																												))
	#### filter low number repeats
	names(sort(table(tmp$repeatclass),decreasing = T))
	list1 <- c("Alu","Simple_repeat","L1","MIR","L2","ERV1","Low_complexity","ERVL-MaLR","hAT-Charlie","ERVL","TcMar-Tigger","CR1","ERVK","significant L1","SVA","Satellite","nonL1PA LINE1","L1PA1-3","other L1PA")
	list <- c(list1)
	tmp1 <- tmp %>% filter(repeatclass %in% list)
	tmp1 <- tmp1 %>% 
		mutate(repeatclass = ifelse(str_detect(repeatclass,"L1"),repeatclass,X6),
					 repeatclass = str_replace(repeatclass,"^LINE$","other LINE"),
					 ) 
	table(tmp1$repeatclass)
	tmp1
	
	### create name with number
	numberdf <- tmp1 %>% group_by(repeatclass) %>% summarize(count=n(),median=median(log2FoldChange)) %>% mutate(group2=paste0(repeatclass,"\n","(",count,")"))
	
	genedf <- gene %>% mutate(group2= paste0("gene","\n","(",nrow(gene),")")) %>% 
		select(group2,log2FoldChange)
	
	# get order by median log2FC
	genenMedian <- genedf %>% group_by(group2) %>% summarize(median = median(log2FoldChange)) 
	tmp2 <- numberdf %>% select(median,group2) %>% rbind(genenMedian)
	levels <- arrange(tmp2,median) %$% group2
	# data frame used for plotting, combine number data frmae
	plotdf <- left_join(tmp1, numberdf %>% dplyr::select(repeatclass,group2))
	# order of plotting

	combinedf <- plotdf %>% select(group2,log2FoldChange) %>% rbind(genedf)
	combinedf$group2  <- factor(combinedf$group2,levels=rev(levels))
	
	
	p2 <- ggplot(combinedf,aes(x=group2,y=log2FoldChange))+
		geom_boxplot(width=0.6,outlier.size = 0,outlier.alpha = 0,notch = T)+
		geom_hline(yintercept = 0,color="red",alpha=0.5,size=0.2)+
		scale_fill_manual(values=c("#4DBBD5FF"))+
		scale_y_continuous(limits = c(-1.7,2.7),expand = c(0,0))+
		theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "right",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
	)+
		guides(fill=guide_legend(title=""))+
		xlab("")
p2	
ggsave(plot = p2,"individual_expression_20221212.pdf",height=2.5,width = 4.5,useDingbats = F)
