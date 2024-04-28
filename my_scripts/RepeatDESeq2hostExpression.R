#!/LiuLab/software/miniconda3/bin/Rscript
library(tidyverse)
library(magrittr)
theme_set(theme_bw())
rm(list=ls())
filter <- dplyr::filter
select <- dplyr::select
desc <- dplyr::desc
reduce <- purrr::reduce
readr::local_edition(1)
options(stringsAsFactors = FALSE)
set.seed(629)
#### three input;
### 1 annotated bed file;
### 2 gene DESeq2 file;
### 3 output file name
args <- commandArgs(trailingOnly = TRUE)
anno_bed <- args[1]
gene_DESeq <- args[2]
output <- args[3]
# 
# anno_bed <- "/analysis/lixf/RNA/human/hESC/KO/SLTM/NL_180929/6-LINE1M/KO/KO_p_0.05.all.annotated.bed"
# 
# gene_DESeq <- "/analysis/lixf/RNA/human/hESC/KO/SLTM/NL_180929/5-DESeq2/KnownGene/Combine_KnownGene.Deseq2.SLTM_KO_vs_SAFE.result.txt"

anno <- read_tsv(anno_bed)
DESeq <- read_tsv(gene_DESeq) 
### same column name
colnames(DESeq)[1] <- "host_gene"

host_gene <- anno$host_gene
### dealing with ";"
all_symbol <- unlist(str_split(host_gene,";")) %>%
	as.data.frame() %>% 
	distinct()
colnames(all_symbol) <- "host_gene"

DESeq_symbol <- left_join(all_symbol,DESeq,by = "host_gene")

list <- list()
n=1
for (gene in host_gene) {
	### test whether the host geen contained ";"
	if (!str_detect(gene,";")) {
		if ( gene %in% DESeq_symbol$host_gene ) {
			list[[n]] <- DESeq_symbol %>% filter(host_gene == gene)
		}
	}
	if (str_detect(gene,";")) {
		gene_split <- str_split_fixed(gene,";",2)
		gene1 <- gene_split[1]
		gene2 <- gene_split[2]
		df1 <- DESeq_symbol %>% filter(host_gene == gene1)
		df2 <- DESeq_symbol %>% filter(host_gene == gene2)
		if(nrow(df1)==0){ 
			df1 <- df1 %>% add_row( host_gene = NA  ) 
			df1[is.na(df1)] <- "-"
			}
		if(nrow(df2)==0){ 
			df2 <- df2 %>% add_row( host_gene = NA  ) 
			df2[is.na(df2)] <- "-"
			}
		df <- t(paste(df1,df2,sep=";")) %>% as.data.frame()
		colnames(df) <- colnames(df1)
		list[[n]] <- df
	}	
	n=n+1
} 

gene_anno_df <- reduce(list,rbind) %>% 
	select(host_gene,log2FoldChange,pvalue,padj) %>% 
	distinct()
colnames(gene_anno_df) <- c("host_gene","gene_log2FoldChange","gene_pvalue","gene_padj"    )
### combine 
total_anno <- left_join(anno,gene_anno_df)
write_tsv(x = total_anno,path = output)
