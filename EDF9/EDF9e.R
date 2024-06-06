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

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/Nature_revision/mouseCRISPRai/Volcano_MAplot/")
file <- c("/analysis/lixf/CRISPRa_i/mESC/RNA/BLY_20230326/CRISPRa/5-DESeq2/RepeatMasker/Combine_Deseq2.merge_L1_all_combine_vs_merge_SAFE.result.txt",
					"/analysis/lixf/CRISPRa_i/mESC/RNA/BLY_20230424_CRISPRi_con/5-DESeq2/RepeatMasker/Combine_Deseq2.merge_L1MdTf_vs_merge_SAFE.result.txt",
					"/analysis2/lixf/mouse_embryo/ShenLab_TaoYB_L1AMO/added_data/5-DESeq2/RepeatMasker/Combine_Deseq2.L1_vs_SCR.result.txt")
file
names(file) <- c("CRISPRa","CRISPRi","L1 AMO")
file
list <- imap(file,~fread(.x) %>% mutate(sample = .y))
df <- reduce(list,rbind)
family <- fread("/LiuLab/reference/Mouse/GRCm38/TE/TE_class/all.name",header =F)
colnames(family) <- c("rowname","class","family")
combine_df <- left_join(df,family)


combine_df <- combine_df %>% mutate(class = ifelse(class %in% c("LINE","SINE","LTR","Retroposon","DNA","SVA","Simple_repeat","Low_complexity"),class, "others")
)

max_padj_log10 <- max(-log10(combine_df$padj[combine_df$padj!=0]))
combine_df <- combine_df %>% 
	mutate(padj = ifelse(padj == 0, 10^-(max_padj_log10 + 50), padj))

LFCcutoff <- 0.1
pvalueL1 <- combine_df %>% filter(family == "L1" & pvalue < 0.1 & abs(log2FoldChange) > LFCcutoff ) %>% group_by(sample) %>%  slice_min(padj, n = 5,with_ties = T) %>% ungroup
#pvalueERVK <- combine_df %>% filter(family == "L1" & padj < 0.01 & abs(log2FoldChange) > LFCcutoff)
#pvalueSVA <- combine_df %>% filter(family == "SVA" & padj < 0.01 )
pvaluedf <- combine_df %>% filter(class %in% c("LTR","DNA","SINE","Retroposon") & pvalue < 0.05 & abs(log2FoldChange) > LFCcutoff)
label_dftmp <- pvaluedf %>% group_by(sample,class) %>% top_bottom(n = 1,wt = log2FoldChange,with_ties = F) %>% ungroup
#label_df <- reduce(list(label_dftmp,pvalueSVA,pvalueL1),rbind) %>% distinct()
label_df <- rbind(pvalueL1)#,label_dftmp)
colorManual <- c(LINE = '#A73030FF',SINE='#EFC000FF',LTR='#0073C2FF',Retroposon='#868686FF',"DNA"='#7AA6DCFF',Simple_repeat = "#003C67FF","Low_complexity" = "#3B3B3BFF",others = "black")
size1 <- 0.8
size2 <- 0.5
sizeManual <- c(LINE = size1,SINE=size1,LTR=size1,Retroposon=size1,"DNA"=size1,Simple_repeat = size2,Low_complexity = size2,others = size2)
alphaManual <- c(LINE = 0.8,SINE=0.8,LTR=0.8,Retroposon=0.8,"DNA"=0.8,Simple_repeat = 0.8,Low_complexity = 0.8,others = 0.2)
combine_df %>% pull(sample) %>% unique
p <- combine_df %>% 
	ungroup %>% 
	#mutate(sample = factor(sample,levels = c("CRISPRa","CRISPRi","AMO") %>% rev)) %>% 
	mutate(sample = fct_inorder(sample)) %>% 
	ggplot(aes(x = log2FoldChange, y = -log10(pvalue),color = class,size = class,alpha = class) ) +
	geom_point(stroke	 = 0)+
	scale_color_manual(values = colorManual)+
	scale_size_manual(values = sizeManual)+
	scale_alpha_manual(values = alphaManual)+
	#	ylim(0,max_padj_log10+50)+
	geom_hline(yintercept = 2,color = "red",alpha = 0.7,linetype = 2,size = 0.15)+
	ggrepel::geom_text_repel(data = label_df,aes(label = rowname),size = 6/ggplot2::.pt,show.legend = FALSE, segment.size = 0.05)+
	facet_wrap(~sample,scales = "free",nrow = 1)+
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				text=element_text(size= 8,family = "ArialMT",color = "black"),
				axis.title = element_text(size = 6,family = "ArialMT",color = "black"),
				plot.subtitle = element_text(hjust = 1),
				axis.text.x = element_text(angle=0,vjust=0),
				plot.margin = margin(0.3,1,0.3,0.3, "cm"),
				strip.placement = "outside",
				strip.background = element_blank(),
				legend.title = element_text(size = 6,family = "ArialMT",color = "black"),
				legend.position = "right",
				legend.key.width= unit(6, 'pt'),
				legend.key.height= unit(4, 'pt'),
	)+guides(color=guide_legend(title="family",override.aes = list(size=1)),alpha="none",size="none",text = "none")
pp <- p+theme(legend.position = "bottom")
ggsave(plot = pp,"./merge_volcanoPlot_bottom_legend_high.pdf",width = 4.0,height  = 2.95,useDingbats = F)
