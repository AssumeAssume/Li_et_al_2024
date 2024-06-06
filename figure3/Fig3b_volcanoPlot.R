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

setwd("/analysis/lixf/CRISPRa_i/NCCIT/RNA//CRISPRi/5-DESeq2/RepeatMasker")
file <- c("/analysis/lixf/CRISPRa_i/NCCIT/RNA//CRISPRi/5-DESeq2/RepeatMasker/Combine_add_RepeatMasker.Deseq2.L1_mix_vs_SAFE.result.txt"
)
names(file) <- c("L1_mix")
list <- imap(file,~fread(.x) %>% mutate(sample = .y))
df <- reduce(list,rbind)
family <- fread("/LiuLab/reference/Human/GRCh38//TE/Repeat_Masker_hg38.class",header =F)
colnames(family) <- c("rowname","class","family")
combine_df <- left_join(df,family)

# selected_df <- combine_df %>% filter(family %in% c("L1","ERVL","ERVK","ERV1","ERVL-MaLR","SVA")) %>% filter(padj< 0.01) %>% filter(abs(log2FoldChange) > 0.5)
# subfamilyName <- selected_df %>% pull(rowname) %>% unique
# 
# longdf <- combine_df %>% filter(rowname %in% subfamilyName)
# longdf <- longdf %>% select(rowname,log2FoldChange,sample) 
# matdf <- dcast(longdf,rowname~sample  ,value.var = "log2FoldChange" )
# mat <- matdf[,-1] %>% as.matrix()
# rownames(mat) <- matdf$rowname
# anno <- matdf %>% left_join(family) %>% select(rowname,family)
# anno_col <- anno %>% select(family)
# rownames(anno_col) <- anno$rowname
# annoCol<-list(family=c(L1 = '#e41a1c',ERV1='#377eb8',ERVL='#4daf4a',ERVK='#984ea3',"ERVL-MaLR"='#ff7f00',SVA = "#ffff33"))
# pdf("./logFC_heatmap.combine.pdf",height =  10,width  = 5)
# pheatmap(mat,
# 				 annotation_row = anno_col,
# 				 color = colorRampPalette(c("blue","white","red"))(300),
# 				 annotation_colors = annoCol,
# 				 cluster_cols = F,
# 				 fontsize = 5,
# 				 cellwidth = 5,cellheight = 5)
# dev.off()

#combine_df <- combine_df %>% filter(sample != "RAB7A")
combine_df <- combine_df %>% mutate(class = ifelse(class %in% c("LINE","SINE","LTR","Retroposon","DNA","SVA","Simple_repeat","Low_complexity"),class, "others")
)

max_padj_log10 <- max(-log10(combine_df$padj[combine_df$padj!=0]))
combine_df <- combine_df %>% 
	mutate(padj = ifelse(padj == 0, 10^-(max_padj_log10 + 50), padj))

LFCcutoff <- 1
pvalueL1 <- combine_df %>% filter(family == "L1" & padj < 0.01 & abs(log2FoldChange) > LFCcutoff ) %>% slice_min(padj, n = 5,with_ties = T)
#pvalueERVK <- combine_df %>% filter(family == "L1" & padj < 0.01 & abs(log2FoldChange) > LFCcutoff)
#pvalueSVA <- combine_df %>% filter(family == "SVA" & padj < 0.01 )
pvaluedf <- combine_df %>% filter(class %in% c("LTR","DNA","SINE","Retroposon") & padj < 0.01 & abs(log2FoldChange) > LFCcutoff)
label_dftmp <- pvaluedf %>% group_by(class) %>% top_bottom(n = 1,wt = log2FoldChange,with_ties = F) %>% ungroup
#label_df <- reduce(list(label_dftmp,pvalueSVA,pvalueL1),rbind) %>% distinct()
label_df <- pvalueL1
colorManual <- c(LINE = '#A73030FF',SINE='#EFC000FF',LTR='#0073C2FF',Retroposon='#868686FF',"DNA"='#7AA6DCFF',Simple_repeat = "#003C67FF","Low_complexity" = "#3B3B3BFF",others = "black")
size1 <- 0.8
size2 <- 0.5
sizeManual <- c(LINE = size1,SINE=size1,LTR=size1,Retroposon=size1,"DNA"=size1,Simple_repeat = size2,Low_complexity = size2,others = size2)
alphaManual <- c(LINE = 0.8,SINE=0.8,LTR=0.8,Retroposon=0.8,"DNA"=0.8,Simple_repeat = 0.8,Low_complexity = 0.8,others = 0.2)

p <- combine_df %>% 
	ggplot(aes(x = log2FoldChange, y = -log10(padj),color = class,size = class,alpha = class) ) +
	geom_point(stroke	 = 0)+
	scale_color_manual(values = colorManual)+
	scale_size_manual(values = sizeManual)+
	scale_alpha_manual(values = alphaManual)+
	ylim(0,max_padj_log10+50)+
	geom_hline(yintercept = 2,color = "red",alpha = 0.7,linetype = 2,size = 0.15)+
	ggrepel::geom_text_repel(data = label_df,aes(label = rowname),size = 6/ggplot2::.pt,show.legend = FALSE, segment.size = 0.05)+
	facet_wrap(~sample,scales = "free",ncol = 4)+
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
p
ggsave(plot = p,"./volcanoPlot.pdf",width = 2.9,height  = 1.7,useDingbats = F)
