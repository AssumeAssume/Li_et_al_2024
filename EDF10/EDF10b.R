library(tidyverse)
library(magrittr)
library(data.table)
library(glue)
library(ggExtra)
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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine//mouseCRISPRai/Gene_cloest_L1//")



geneList <- fread("/analysis/lixf/CRISPRa_i/mESC/RNA//public_gene_sets/2c.txt")
colnames(geneList)[1] <- "GeneName"
DBTMEE <- fread("/analysis/lixf/CRISPRa_i/mESC/RNA//public_gene_sets/DBTMEE/tables/cluster_gene_v2.tsv")
DBTMEE <- DBTMEE %>% rename(geneset = Cluster,GeneName = Gene) %>% select(geneset,GeneName)
majorZGA <- DBTMEE %>% filter(geneset=="Major ZGA")
minorZGA <- DBTMEE %>% filter(geneset=="Minor ZGA")

DistanceDf <- fread("../Gene_6K_L1.distance.txt") 
addedGene <- DistanceDf %>% filter(str_detect(V4,"Zscan4|MERVL|MT2_|BB287469|BC080695|Gm2076|Gm21304|Gm21312|Gm4027|Gm4340|Gm5039|Gm7978|Gm8300|Gm8764|Gm8994|Otud6a|Pramef25|Sp110|Usp17|Obox4-ps12|Obox4-ps16|Obox4-ps18|Obox4-ps23|Obox4-ps35")) %>% select(V4) %>% rename(GeneName = V4)
geneList <- rbind(geneList,addedGene)


### 2c gene sets  ----------------------------------------------------------------------------------------------------------------------------------

# crispra -----------------------------------------------------------------------------------------------------------------------------

expA <- fread("/analysis/lixf/CRISPRa_i/mESC/RNA/combine_batch/5-DESeq2/KnownGene/Combine_Deseq2.combine_CRISPRa_L1_mix_vs_combine_CRISPRa_SAFE.result.txt")
colnames(expA)[1] <- "V4"
DistanceDf %>% 
	filter(V13 != -1) %>% 
	mutate(group = ifelse(V4 %in% geneList$GeneName,"2C","Non-2C")) %>% 
	ggplot(aes(x = group,y = log2(abs(V13+1))))+
	geom_boxplot()+
	stat_compare_means()
plotdf <- DistanceDf %>% 
	filter(V13 != -1) %>% 
	mutate(group = ifelse(V4 %in% geneList$GeneName,"2C-specific genes","others")) %>% 
	full_join(expA) %>% 
	filter(! (is.na(V13)|is.na(log2FoldChange)) ) 

label <- plotdf %>% filter( V4 %in% geneList$GeneName ) %>%
	mutate(labelName = case_when(str_detect(V4,"Usp17l") ~ "Usp17l",
															 str_detect(V4,"Zscan4") & !str_detect(V4,"Zscan4-ps") ~ "Zscan4",
															 str_detect(V4,"Obox4-ps") ~ "Obox4-ps",
															 str_detect(V4,"Gm8300|Gm4027|Gm16381|Gm21319|BB287469|Gm2022|Gm5662|Gm5039") ~ "Eif1a-like",
															 str_detect(V4,"Tdpoz|Gm9125|Gm10697|Gm4858") ~ "Tdpoz",
															 str_detect(V4,"Pramef25|BC080695|Gm13057|Gm13040|Gm13043|Gm3106|Gm3139|Gm16513|D5Ertd577e|Gm3286|E330014E10Rik|Gm3259|Gm7978|Gm7982|Gm6367|BC080696|AA792892|A430089I19Rik|BC061212|C87414") ~ "Prame-like",
															 str_detect(V4,"Gm20767|Tcstv|AF067061|BC147527|AF067063") ~ "Tcstv",
															 str_detect(V4,"Gm4340|Gm21304|GM20765|Gm21312|Gm8764|Gm6763|Gm21293") ~ "chr10 Unknown Cluster",
															 TRUE ~ "else"
															 )) %>% 
	filter(labelName != "else") %>% 
	group_by(labelName) %>% 
	filter( log2FoldChange == max(log2FoldChange)) %>% 
	ungroup
#fwrite(label,"label_tmp.tsv",sep = '\t')

#label %>% pull(V4)
max=1.8e-06
plotdf %>% arrange(desc(group)) %>% 
	ggplot(aes(x = V13,color = group,group=group,alpha=group,size=group))+
	geom_density()+
	scale_alpha_manual(values = c("2C-specific genes"=1,"others"=0.5))+
	scale_size_manual(values = c("2C-specific genes"=0.6,"others"=0.3))+
	scale_color_manual(values = c("2C-specific genes"="#A83030","others"="grey"))+
	scale_y_continuous(breaks= max,labels = max)
ggsave("CRISPRa_density.pdf",width = 5,height = 1)
p1 <- plotdf %>% arrange(desc(group)) %>% 
	ggplot(aes(x = V13, y = log2FoldChange,color = group,group=group,alpha=group,size=group))+
	geom_point(stroke= 0.1)+
	scale_alpha_manual(values = c("2C-specific genes"=1,"others"=0.5))+
	scale_size_manual(values = c("2C-specific genes"=0.6,"others"=0.3))+
	scale_color_manual(values = c("2C-specific genes"="#A83030","others"="grey"))+
	ylim ( c(-1,2.5))+
	scale_x_continuous(breaks =seq(-10000000,10000000,1000000),limits = c(-5000000,5000000),
										 label = seq(-10000000,10000000,1000000)/1000000)+
	ggrepel::geom_text_repel(data=label,aes(label = V4),size = ggplot2::.pt,min.segment.length = 0,segment.size = 0.1,
													 nudge_x      = -4000000,
													 direction    = "y",
													 hjust        = "left")+
	geom_vline(xintercept = -512000,color = "black",alpha = 0.7,linetype = 2)+
	geom_vline(xintercept = 512000,color = "black",alpha = 0.7,linetype = 2)+
	geom_hline(yintercept = 0,color = "black",alpha = 0.7,linetype = 2)+
	xlab("TSS to nearest FL-L1")+ylab("LFC CRISPRa")+
	theme(legend.position = "bottom")+theme(panel.grid.major = element_blank(),
																					panel.grid.minor = element_blank(),
																					text=element_text(size= 8,family = "ArialMT",color = "black"),
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
	)
#ph <- ggMarginal(p1, type = "boxplot",margins = "both",groupColour  = T,groupFill = T)
pd <- ggMarginal(p1, type = "density",margins = "x",groupColour  = T,groupFill = T,width = 0.4)
#ph$grobs[ph$layout$name == "topMargPlot"] <- pd$grobs[ph$layout$name == "topMargPlot"]
#ph
ggsave(plot=pd,"CRISPRa.pdf",width = 2.5,height = 3,useDingbats= F)

#ggsave(plot=pd,"CRISPRa_large.pdf",width = 5.5,height = 5,useDingbats= F)

plotdf %>% mutate(distance = abs(V13)) %>% wilcox_test(distance~group)

plotdf %>% arrange(desc(group)) %>% 
	ggplot(aes(x = group, y = log2FoldChange,color = group))+
	geom_boxplot(width = 0.5,outlier.alpha = 0)+
	scale_color_manual(values = c("2C-specific genes"="#A83030","others"="grey"))+
	theme_classic()+
	geom_hline(yintercept = 0,color = "black",alpha = 0.7,linetype = 2)+
	ylim(c(-0.5,1.5))+
	xlab("")+
	scale_x_discrete(label = c("",""))+
	ggpubr::stat_compare_means(label.y =1.0,size= 6/ggplot2::.pt, aes(label = paste0("p = ", after_stat(p.format))))+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "none",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt"),
				axis.line =  element_line(colour = "black", size = 0.25)
	)
ggsave("CRISPRa_LFC_boxplot.pdf",width = 1.5,height = 2)

# CRISPRi -----------------------------------------------------------------------------------------------------------------------------


expI <- fread("/analysis/lixf/CRISPRa_i/mESC/RNA/5-DESeq2/KnownGene/Combine_Deseq2.merge_L1_all_combine_vs_merge_SAFE.result.txt")
colnames(expI)[1] <- "V4"
plotdf <- DistanceDf %>% 
	filter(V13 != -1) %>% 
	mutate(group = ifelse(V4 %in% geneList$GeneName,"2C-specific genes","others")) %>% 
	full_join(expI) %>% 
	filter(! (is.na(V13)|is.na(log2FoldChange)) ) 
label <- plotdf %>% filter( V4 %in% geneList$GeneName ) %>%
	mutate(labelName = case_when(str_detect(V4,"Usp17l") ~ "Usp17l",
															 str_detect(V4,"Zscan4")& !str_detect(V4,"Zscan4-ps") ~ "Zscan4",
															 str_detect(V4,"Obox4-ps") ~ "Obox4-ps",
															 str_detect(V4,"Gm8300|Gm4027|Gm16381|Gm21319|BB287469|Gm2022|Gm5662|Gm5039") ~ "Eif1a-like",
															 str_detect(V4,"Tdpoz|Gm9125|Gm10697|Gm4858") ~ "Tdpoz",
															 str_detect(V4,"Pramef25|BC080695|Gm13057|Gm13040|Gm13043|Gm3106|Gm3139|Gm16513|D5Ertd577e|Gm3286|E330014E10Rik|Gm3259|Gm7978|Gm7982|Gm6367|BC080696|AA792892|A430089I19Rik|BC061212|C87414") ~ "Prame-like",
															 str_detect(V4,"Gm20767|Tcstv|AF067061|BC147527|AF067063") ~ "Tcstv",
															 str_detect(V4,"Gm4340|Gm21304|GM20765|Gm21312|Gm8764|Gm6763|Gm21293") ~ "chr10 Unknown Cluster",
															 TRUE ~ "else"
	)) %>% 
	filter(labelName != "else") %>% 
	group_by(labelName) %>% 
	filter( log2FoldChange == min(log2FoldChange)) %>% 
	ungroup
max=1.7e-06
plotdf %>% arrange(desc(group)) %>% 
	ggplot(aes(x = V13,color = group,group=group,alpha=group,size=group))+
	geom_density()+
	scale_alpha_manual(values = c("2C-specific genes"=1,"others"=0.5))+
	scale_size_manual(values = c("2C-specific genes"=0.6,"others"=0.3))+
	scale_color_manual(values = c("2C-specific genes"="#A83030","others"="grey"))+
	scale_y_continuous(breaks= max,labels = max)
ggsave("CRISPRi_density.pdf",width = 5,height = 1)

p2 <- plotdf %>% arrange(desc(group)) %>% 
	ggplot(aes(x = V13, y = log2FoldChange,color = group,group=group,alpha=group,size=group))+
	geom_point(stroke= 0.1)+
	scale_alpha_manual(values = c("2C-specific genes"=1,"others"=0.5))+
	scale_size_manual(values = c("2C-specific genes"=0.6,"others"=0.3))+
	scale_color_manual(values = c("2C-specific genes"="#A83030","others"="grey"))+
	ylim ( c(-1,0.6))+
	scale_x_continuous(breaks =seq(-10000000,10000000,1000000),limits = c(-5000000,5000000),
										 label = seq(-10000000,10000000,1000000)/1000000)+
	ggrepel::geom_text_repel(data=label,aes(label = V4),size = ggplot2::.pt,min.segment.length = 0,segment.size = 0.1,
													 nudge_x      = -4000000,
													 direction    = "y",
													 hjust        = "left")+
	geom_vline(xintercept = -512000,color = "black",alpha = 0.7,linetype = 2)+
	geom_vline(xintercept = 512000,color = "black",alpha = 0.7,linetype = 2)+
	geom_hline(yintercept = 0,color = "black",alpha = 0.7,linetype = 2)+
	xlab("TSS to nearest FL-L1")+ylab("LFC CRISPRa")+
	theme(legend.position = "bottom")+theme(panel.grid.major = element_blank(),
																					panel.grid.minor = element_blank(),
																					text=element_text(size= 8,family = "ArialMT",color = "black"),
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
	)
#ph <- ggMarginal(p1, type = "boxplot",margins = "both",groupColour  = T,groupFill = T)
pd <- ggMarginal(p2, type = "density",margins = "x",groupColour  = T,groupFill = T,width = 0.4)
#ph$grobs[ph$layout$name == "topMargPlot"] <- pd$grobs[ph$layout$name == "topMargPlot"]
#ph
ggsave(plot=pd,"CRISPRi.pdf",width = 2.5,height = 3,useDingbats= F)
plotdf %>% mutate(distance = abs(V13)) %>% wilcox_test(distance~group)

plotdf %>% arrange(desc(group)) %>% 
	ggplot(aes(x = group, y = log2FoldChange,color = group))+
	geom_boxplot(width = 0.5,outlier.alpha = 0)+
	scale_color_manual(values = c("2C-specific genes"="#A83030","others"="grey"))+
	theme_classic()+
	geom_hline(yintercept = 0,color = "black",alpha = 0.7,linetype = 2)+
	ylim(c(-0.4,0.3))+
	xlab("")+
	scale_x_discrete(label = c("",""))+
	ggpubr::stat_compare_means(label.y =-.3,size= 6/ggplot2::.pt, aes(label = paste0("p = ", after_stat(p.format))))+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "none",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt"),
				axis.line =  element_line(colour = "black", size = 0.25)
	)
ggsave("CRISPRi_LFC_boxplot.pdf",width = 1.5,height = 2)




# mESC ASO ----------------------------------------------------------------------------------------------------------------------------


expASO <- fread("/analysis2/lixf/mouse_embryo/added_data/5-DESeq2/KnownGene/Combine_Deseq2.L1_vs_SCR.result.txt")
colnames(expASO)[1] <- "V4"
plotdf <- DistanceDf %>% 
	filter(V13 != -1) %>% 
	mutate(group = ifelse(V4 %in% geneList$GeneName,"2C-specific genes","others")) %>% 
	full_join(expASO) %>% 
	filter(! (is.na(V13)|is.na(log2FoldChange)) ) 
max=1.25e-06
plotdf %>% arrange(desc(group)) %>% 
	ggplot(aes(x = V13,color = group,group=group,alpha=group,size=group))+
	geom_density()+
	scale_alpha_manual(values = c("2C-specific genes"=1,"others"=0.5))+
	scale_size_manual(values = c("2C-specific genes"=0.6,"others"=0.3))+
	scale_color_manual(values = c("2C-specific genes"="#A83030","others"="grey"))+
	scale_y_continuous(breaks= max,labels = max)
ggsave("mESC_ASO_density.pdf",width = 5,height = 1)


label <- plotdf %>% filter( V4 %in% geneList$GeneName ) %>%
	mutate(labelName = case_when(str_detect(V4,"Usp17l") ~ "Usp17l",
															 str_detect(V4,"Zscan4")& !str_detect(V4,"Zscan4-ps") ~ "Zscan4",
															 str_detect(V4,"Obox4-ps") ~ "Obox4-ps",
															 str_detect(V4,"Gm8300|Gm4027|Gm16381|Gm21319|BB287469|Gm2022|Gm5662|Gm5039") ~ "Eif1a-like",
															 str_detect(V4,"Tdpoz|Gm9125|Gm10697|Gm4858") ~ "Tdpoz",
															 str_detect(V4,"Pramef25|BC080695|Gm13057|Gm13040|Gm13043|Gm3106|Gm3139|Gm16513|D5Ertd577e|Gm3286|E330014E10Rik|Gm3259|Gm7978|Gm7982|Gm6367|BC080696|AA792892|A430089I19Rik|BC061212|C87414") ~ "Prame-like",
															 str_detect(V4,"Gm20767|Tcstv|AF067061|BC147527|AF067063") ~ "Tcstv",
															 str_detect(V4,"Gm4340|Gm21304|GM20765|Gm21312|Gm8764|Gm6763|Gm21293") ~ "chr10 Unknown Cluster",
															 TRUE ~ "else"
	)) %>% 
	filter(labelName != "else") %>% 
	group_by(labelName) %>% 
	filter( log2FoldChange == min(log2FoldChange)) %>% 
	ungroup
p3 <- plotdf %>% arrange(desc(group)) %>% 
	ggplot(aes(x = V13, y = log2FoldChange,color = group,group=group,alpha=group,size=group))+
	geom_point(stroke= 0.1)+
	scale_alpha_manual(values = c("2C-specific genes"=1,"others"=0.5))+
	scale_size_manual(values = c("2C-specific genes"=0.6,"others"=0.3))+
	scale_color_manual(values = c("2C-specific genes"="#A83030","others"="grey"))+
	ylim ( c(-1.5,1))+
	scale_x_continuous(breaks =seq(-10000000,10000000,1000000),limits = c(-5000000,5000000),
										 label = seq(-10000000,10000000,1000000)/1000000)+
	ggrepel::geom_text_repel(data=label,aes(label = V4),size = ggplot2::.pt,min.segment.length = 0,segment.size = 0.1,
													 nudge_x      = -4000000,
													 direction    = "y",
													 hjust        = "left")+
	geom_vline(xintercept = -512000,color = "black",alpha = 0.7,linetype = 2)+
	geom_vline(xintercept = 512000,color = "black",alpha = 0.7,linetype = 2)+
	geom_hline(yintercept = 0,color = "black",alpha = 0.7,linetype = 2)+
	xlab("TSS to nearest FL-L1")+ylab("LFC mESC ASO")+
	theme(legend.position = "bottom")+theme(panel.grid.major = element_blank(),
																					panel.grid.minor = element_blank(),
																					text=element_text(size= 8,family = "ArialMT",color = "black"),
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
	)
p3
#ph <- ggMarginal(p1, type = "boxplot",margins = "both",groupColour  = T,groupFill = T)
pd <- ggMarginal(p3, type = "density",margins = "x",groupColour  = T,groupFill = T,width = 0.4)

#ph$grobs[ph$layout$name == "topMargPlot"] <- pd$grobs[ph$layout$name == "topMargPlot"]
#ph
ggsave(plot=pd,"ASO_mESC.pdf",width = 2.5,height = 3,useDingbats= F)
plotdf %>% wilcox_test(log2FoldChange~group)
plotdf %>% mutate(distance = abs(V13)) %>% wilcox_test(distance~group)
plotdf %>% arrange(desc(group)) %>% 
	ggplot(aes(x = group, y = log2FoldChange,color = group))+
	geom_boxplot(width = 0.5,outlier.alpha = 0)+
	scale_color_manual(values = c("2C-specific genes"="#A83030","others"="grey"))+
	theme_classic()+
	geom_hline(yintercept = 0,color = "black",alpha = 0.7,linetype = 2)+
	ylim(c(-0.8,0.6))+
	xlab("")+
	scale_x_discrete(label = c("",""))+
	ggpubr::stat_compare_means(label.y =-.7,size= 6/ggplot2::.pt, aes(label = paste0("p = ", after_stat(p.format))))+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "none",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt"),
				axis.line =  element_line(colour = "black", size = 0.25)
	)
ggsave("mESC_ASO_LFC_boxplot.pdf",width = 1.5,height = 2)



# 2c embryo ---------------------------------------------------------------------------------------------------------------------------


# 2c embryo ---------------------------------------------------------------------------------------------------------------------------




ASO2c <- fread("/analysis/lixf/CRISPRa_i/mESC/RNA/late2C_AMO/5-DESeq2/KnownGene/Combine_removeL1_1Deseq2.L1_vs_SCR.result.txt")
colnames(ASO2c)[1] <- "V4"
plotdf <- DistanceDf %>% 
	filter(V13 != -1) %>% 
	mutate(group = ifelse(V4 %in% geneList$GeneName,"2C-specific genes","others")) %>% 
	full_join(ASO2c) %>% 
	filter(! (is.na(V13)|is.na(log2FoldChange)) ) 

max=2.5e-06
plotdf %>% arrange(desc(group)) %>% 
	ggplot(aes(x = V13,color = group,group=group,alpha=group,size=group))+
	geom_density()+
	scale_alpha_manual(values = c("2C-specific genes"=1,"others"=0.5))+
	scale_size_manual(values = c("2C-specific genes"=0.6,"others"=0.3))+
	scale_color_manual(values = c("2C-specific genes"="#A83030","others"="grey"))+
	scale_y_continuous(breaks= max,labels = max)
ggsave("own_Late2C_0915_ASO_density.pdf",width = 5,height = 1)

label <- plotdf %>% filter( V4 %in% geneList$GeneName ) %>%
	mutate(labelName = case_when(str_detect(V4,"Usp17l") ~ "Usp17l",
															 str_detect(V4,"Zscan4") & !str_detect(V4,"Zscan4-ps")~ "Zscan4",
															 str_detect(V4,"Obox4-ps") ~ "Obox4-ps",
															 str_detect(V4,"Gm8300|Gm4027|Gm16381|Gm21319|BB287469|Gm2022|Gm5662|Gm5039") ~ "Eif1a-like",
															 str_detect(V4,"Tdpoz|Gm9125|Gm10697|Gm4858") ~ "Tdpoz",
															 str_detect(V4,"Pramef25|BC080695|Gm13057|Gm13040|Gm13043|Gm3106|Gm3139|Gm16513|D5Ertd577e|Gm3286|E330014E10Rik|Gm3259|Gm7978|Gm7982|Gm6367|BC080696|AA792892|A430089I19Rik|BC061212|C87414") ~ "Prame-like",
															 str_detect(V4,"Gm20767|Tcstv|AF067061|BC147527|AF067063") ~ "Tcstv",
															 str_detect(V4,"Gm4340|Gm21304|GM20765|Gm21312|Gm8764|Gm6763|Gm21293") ~ "chr10 Unknown Cluster",
															 TRUE ~ "else"
	)) %>% 
	filter(labelName != "else") %>% 
	group_by(labelName) %>% 
	filter( log2FoldChange == min(log2FoldChange)) %>% 
	ungroup
p4 <- plotdf %>% arrange(desc(group)) %>% 
	ggplot(aes(x = V13, y = log2FoldChange,color = group,group=group,alpha=group,size=group))+
	geom_point(stroke= 0.1)+
	scale_alpha_manual(values = c("2C-specific genes"=1,"others"=0.5))+
	scale_size_manual(values = c("2C-specific genes"=0.6,"others"=0.3))+
	scale_color_manual(values = c("2C-specific genes"="#A83030","others"="grey"))+
	ylim ( c(-3,2))+
	scale_x_continuous(breaks =seq(-10000000,10000000,1000000),limits = c(-5000000,5000000),
										 label = seq(-10000000,10000000,1000000)/1000000)+
	ggrepel::geom_text_repel(data=label,aes(label = V4),size = ggplot2::.pt,min.segment.length = 0,segment.size = 0.1,
													 nudge_x      = -4000000,
													 direction    = "y",
													 hjust        = "left")+
	geom_vline(xintercept = -512000,color = "black",alpha = 0.7,linetype = 2)+
	geom_vline(xintercept = 512000,color = "black",alpha = 0.7,linetype = 2)+
	geom_hline(yintercept = 0,color = "black",alpha = 0.7,linetype = 2)+
	xlab("TSS to nearest FL-L1")+ylab("LFC 2C ASO")+
	theme(legend.position = "bottom")+theme(panel.grid.major = element_blank(),
																					panel.grid.minor = element_blank(),
																					text=element_text(size= 8,family = "ArialMT",color = "black"),
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
																					#axis.line =  element_line(colour = "black", size = 0.25)
	)
p4
#ph <- ggMarginal(p1, type = "boxplot",margins = "both",groupColour  = T,groupFill = T)
pd <- ggMarginal(p4, type = "density",margins = "x",groupColour  = T,groupFill = T,width = 0.4)
#ph$grobs[ph$layout$name == "topMargPlot"] <- pd$grobs[ph$layout$name == "topMargPlot"]
pd
plotdf %>% wilcox_test(log2FoldChange~group)
plotdf %>% mutate(distance = abs(V13)) %>% wilcox_test(distance~group)
ggsave(plot=pd,"own_Late2C_0915_ASO.pdf",width = 2.5,height = 3,useDingbats= F)

plotdf %>% arrange(desc(group)) %>% 
	ggplot(aes(x = group, y = log2FoldChange,color = group))+
	geom_boxplot(width = 0.5,outlier.alpha = 0)+
	scale_color_manual(values = c("2C-specific genes"="#A83030","others"="grey"))+
	theme_classic()+
	geom_hline(yintercept = 0,color = "black",alpha = 0.7,linetype = 2)+
	ylim(c(-2.2,1.5))+
	xlab("")+
	scale_x_discrete(label = c("",""))+
	ggpubr::stat_compare_means(label.y =-.7,size= 6/ggplot2::.pt, aes(label = paste0("p = ", after_stat(p.format))))+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.position = "none",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
				axis.ticks = element_line(colour = "black", size = 0.25),
				axis.ticks.length=unit(1.5, "pt"),
				axis.line =  element_line(colour = "black", size = 0.25)
	)
ggsave("own_late2C_0915_ASO_LFC_boxplot.pdf",width = 1.5,height = 2)


# combine -----------------------------------------------------------------------------------------------------------------------------

((p1|p2|p3|p4) + plot_layout(guides =  'collect') )&
	theme(legend.position='bottom')
ggsave("Combined_point_plot.pdf",width = 12,height = 3.3)
