library(tidyverse)
library(magrittr)
library(scales)
theme_set(theme_bw())
rm(list=ls())
filter <- dplyr::filter
select <- dplyr::select
reduce <- purrr::reduce
options(stringsAsFactors = FALSE)
set.seed(629)

read_mutate <- function(DESeq2RepeatMasker,prefix){
  tmp <- read_tsv(file = DESeq2RepeatMasker) %>% select(rowname,log2FoldChange,padj,pvalue)
  tmp$padj[is.na(tmp$padj)] <- 1 
  tmp %>% mutate(sample=rep(prefix,length(rowname)))
}

j="RepeatMasker"
#genelist=c("ZNF638","DBR1","DIS3","LSM11","CTBP1","NDUFA1","YTHDC1","THOC1","NDUFB8","YY1")


LINE1_factor <- read.table("/LiuLab/reference/Mouse/GRCm38/TE/LINE_1_factor.txt",col.names = F) %>% as.list() %>% unlist


dir.create("/analysis2/lixf/m6A/mouse/RNA/m6A_combine")
setwd("/analysis2/lixf/m6A/mouse/RNA/m6A_combine")

# df <- reduce(list,rbind)
# fwrite(df,"combine_RepeatMasker_DESeq2_df.tsv",sep = '\t')
df <- fread("./combine_RepeatMasker_DESeq2_df.tsv")
df$sample <- factor(df$sample)
levels(df$sample)
select <- c("Chelmicki et al. 2021 Mettl14" ,"Chelmicki et al. 2021 Mettl3" ,"Chen et al. 2021 Mettl3","Xu et al. 2021 Mettl3")
L1_select <- c("L1Md|L1_Mus|L1_Mur|L1VL")
plotdf <- df %>%
  filter(sample %in% select) %>% 
  filter( str_detect(rowname,"^L1")) %>% 
  filter(str_detect(rowname,L1_select)) %>% 
  mutate(LINE1=factor(rowname,levels=LINE1_factor)) %>% 
  mutate(pvalue_group = case_when(pvalue <= 0.01 ~ "<= 0.01",
                                  pvalue > 0.01 & pvalue <= 0.05  ~ "<= 0.05",
                                  TRUE ~ "> 0.05")
  ) %>% 
  arrange(desc(log2FoldChange)) 

sizeManual <- c("<= 0.01" = 3, "<= 0.05" = 2.3 , "> 0.05" = 1.0)
max <- plotdf$log2FoldChange %>% max()
min <- plotdf$log2FoldChange %>% min()
range <- max - min
midPoint <- 1 - max/range


#### LINE-1 ------
p <- plotdf %>%
  ggplot(aes(x=LINE1, y=sample, size=pvalue_group, fill=log2FoldChange)) +
  geom_point(alpha=1, shape=21, color="black",stroke=0.005) +
  scale_size_manual(values = sizeManual ,name="pvalue") +
  scale_fill_gradientn( colours = c("darkblue","#5A4AA5","white","#A87A70",muted("red")),
                        values = c(0, midPoint - (1/range) ,midPoint,midPoint + (1/range) ,1),)+
  # theme_minimal() +
  theme(panel.grid.major = element_blank(),
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
        axis.ticks.length=unit(1.5, "pt")
  )+
  guides(fill = guide_colorbar(title.position = "top",
                               title.hjust = .5,
                               barwidth = unit(5, "lines"),
                               barheight = unit(.5, "lines")),
         size = guide_legend(title.position = "top",
                             title.hjust = .5,
                             #barwidth = unit(10, "lines"),
                             #barheight = unit(.5, "lines")
         )
  )+
  xlab("")+ylab("") 
ggsave(plot = p,"./selected_pvalue_LINE1_foldchange.pdf",width = 3.5,height =2,useDingbats=FALSE)  
#guides(fill=guide_legend(title="padj"))+
p2 <- p+coord_flip()
ggsave(plot=p2,"./flip_selected_pvalue_LINE1_foldchange.pdf",width = 2.0,height =3.8,useDingbats=FALSE)  

