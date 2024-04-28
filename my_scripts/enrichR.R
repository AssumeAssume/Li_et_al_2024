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
options(stringsAsFactors = FALSE)
set.seed(629)
source("/analysis/lixf/my_script/my_R_function.R")

library(enrichR)
setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
head(dbs)
dbs %>% filter(str_detect(libraryName,"GO"))

dbs_GO <- dbs %>% filter(str_detect(libraryName,"GO")) %>% filter(str_detect(libraryName,"2021")) %>% pull(libraryName)

dbs_domains <- dbs %>% filter(str_detect(libraryName,"Domains")) %>% pull(libraryName)

dbs_celltype <- dbs %>% filter(libraryName %in% c("PanglaoDB_Augmented_2021","ESCAPE","ProteomicsDB_2020","CellMarker_Augmented_2021","ARCHS4_Cell-lines","CCLE_Proteomics_2020","ARCHS4_Tissues")) %>% pull(libraryName)

dbs_tf <- dbs %>% filter(libraryName %in% c("ChEA_2016","ENCODE_and_ChEA_Consensus_TFs_from","TF_Perturbations_Followed_by_Expression","TRRUST_Transcription_Factors_2019","ENCODE_TF_ChIP-seq_2015","TF-LOF_Expression_from_GEO","ENCODE_Histone_Modifications_2015","Transcription_Factor_PPIs","Epigenomics_Roadmap_HM_ChIP-seq","TRANSFAC_and_JASPAR_PWMs") ) %>% pull(libraryName)

dbs_diseases <- dbs %>% filter(libraryName %in% c("PheWeb_2019","ClinVar_2019","Orphanet_Augmented_2021","PhenGenI_Association_2021","GWAS_Catalog_2019
","UK_Biobank_GWAS_v1","DisGeNET","DSigDB","OMIM_Disease","OMIM_Expanded","VirusMINT","MSigDB_Oncogenic_Signatures","Rare_Diseases_GeneRIF_ARCHS4_Predictions","Rare_Diseases_GeneRIF_Gene_Lists","Rare_Diseases_AutoRIF_Gene_Lists","Achilles_fitness_decrease","dbGaP","") ) %>% pull(libraryName)
dbs_pathway <- dbs %>% filter(libraryName %in% c("BioPlanet_2019","MSigDB_Hallmark_2020","WikiPathway_2021_Human","KEGG_2021_Human","Elsevier_Pathway_Collection","BioCarta_2016","Reactome_2022","CORUM") ) %>% pull(libraryName)

candidate_dbs <- c(dbs_pathway,dbs_GO,dbs_domains,dbs_celltype,dbs_tf,dbs_diseases)



######### 
# setwd("/analysis/lixf/RNA/human/K562/KnockDown/DBR1/WY_20201230/")
# dir.create("./enrichR_enrichment")
# setwd("./enrichR_enrichment/")
DBR1 <- fread("/analysis/lixf/RNA/human/K562/KnockDown/DBR1/WY_20201230/5-DESeq2/KnownGene/Combine_KnownGene.Deseq2.DBR1_vs_SAFE.sig.result.txt")

logfccutoff <- 0.5
padjcutoff <- 0.01

upGene <- DBR1 %>% filter(log2FoldChange > logfccutoff & padj < padjcutoff) %>% pull(rowname)

enriched <- enrichr(upGene, candidate_dbs)

filter_list <- imap(enriched,function(x,name) {
  x %>% filter(Adjusted.P.value < 0.1) %>% mutate(database = rep(name,length(Adjusted.P.value)))
} )
adjDF <- reduce(filter_list, rbind)
GO_pathway_adjDF <- adjDF %>% filter(database %in% c(dbs_GO,dbs_pathway))
top <- GO_pathway_adjDF %>% 
  group_by(database) %>% 
  top_n(-5,Adjusted.P.value)


logfccutoff <- 1
padjcutoff <- 0.01

downGene <- DBR1 %>% filter(log2FoldChange < -logfccutoff & padj < padjcutoff) %>% pull(rowname)
enriched_down <- enrichr(downGene, candidate_dbs)

filter_list.2 <- imap(enriched_down,function(x,name) {
  x %>% filter(Adjusted.P.value < 0.1) %>% mutate(database = rep(name,length(Adjusted.P.value)))
} )
adjDF.2 <- reduce(filter_list.2, rbind)
GO_pathway_adjDF.2 <- adjDF.2 %>% filter(database %in% c(dbs_GO,dbs_pathway))
top.2 <- GO_pathway_adjDF.2 %>% 
  group_by(database) %>% 
  top_n(-5,Adjusted.P.value)
top.2 %>% 
  ggplot(aes(x=Term,y=-log10(Adjusted.P.value)))+
  geom_bar(stat="identity",fill = "#49598a",width = 0.5)+
  coord_flip()+xlab("")+ylab(expression(-log["10"]~"Adj. P value"))+
  lemon::facet_rep_wrap(~database,ncol = 1,scales = "free_y")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text( size=rel(1)),
        axis.title.y = element_text( size=rel(1)),
        plot.title = element_text(hjust = 0.5,size=rel(1),face = "bold"),
        plot.subtitle = element_text(hjust = 1,size=rel(1)),
        axis.text.x = element_text(angle=0,vjust=0),
        axis.text= element_text(size=rel(1)),
        legend.text=element_text(size=rel(1)),
        legend.title = element_text(size=rel(1)),
        plot.margin = margin(0.3,1,0.3,0.3, "cm"),
        strip.text.x = element_text(size = rel(1)),
        strip.placement = "outside",strip.background = element_blank(),
        legend.position = "bottom"
  )+guides(fill=guide_legend(title=""))



top.2 %>% 
  mutate(Term=fct_reorder(Term, -log10(P.value),mean)) %>% 
  ggplot(aes(x=Term,y=-log10(P.value)))+
  geom_bar(stat="identity",fill = "#49598a",width = 0.5)+
  coord_flip()+xlab("")+ylab(expression(-log["10"]~"P.value"))+
  lemon::facet_rep_wrap(~database,ncol = 2,scales = "free_y",repeat.tick.labels = T)+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text( size=rel(1)),
        axis.title.y = element_text( size=rel(1)),
        plot.title = element_text(hjust = 0.5,size=rel(1),face = "bold"),
        plot.subtitle = element_text(hjust = 1,size=rel(1)),
        axis.text.x = element_text(angle=0,vjust=0),
        axis.text= element_text(size=rel(1)),
        legend.text=element_text(size=rel(1)),
        legend.title = element_text(size=rel(1)),
        plot.margin = margin(0.3,1,0.3,0.3, "cm"),
        strip.text.x = element_text(size = rel(1)),
        strip.placement = "outside",strip.background = element_blank(),
        legend.position = "bottom"
  )+guides(fill=guide_legend(title=""))
#ggsave("./enrichR_DBR1_downgene.pdf",width = 20,height = 12)

top %>% 
  mutate(Term=fct_reorder(Term, -log10(P.value),mean)) %>% 
  ggplot(aes(x=Term,y=-log10(P.value)))+
  geom_bar(stat="identity",fill = "#49598a",width = 0.5)+
  coord_flip()+xlab("")+ylab(expression(-log["10"]~"P.value"))+
  lemon::facet_rep_wrap(~database,ncol = 2,scales = "free_y",repeat.tick.labels = T)+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text( size=rel(1)),
        axis.title.y = element_text( size=rel(1)),
        plot.title = element_text(hjust = 0.5,size=rel(1),face = "bold"),
        plot.subtitle = element_text(hjust = 1,size=rel(1)),
        axis.text.x = element_text(angle=0,vjust=0),
        axis.text= element_text(size=rel(1)),
        legend.text=element_text(size=rel(1)),
        legend.title = element_text(size=rel(1)),
        plot.margin = margin(0.3,1,0.3,0.3, "cm"),
        strip.text.x = element_text(size = rel(1)),
        strip.placement = "outside",strip.background = element_blank(),
        legend.position = "bottom"
  )+guides(fill=guide_legend(title=""))
#ggsave("./enrichR_DBR1_upgene.pdf",width = 20,height = 12)
