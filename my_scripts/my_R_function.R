filter <- dplyr::filter
select <- dplyr::select
desc <- dplyr::desc
reduce <- purrr::reduce
require(ggplot2)
theme_set(theme_bw())

############# screening results ##############
# 
L1_sup <- vroom::vroom("/analysis/lixf/sgRNA/human/LINE_screening_combine/suppressor.csv")
L1_act <- vroom::vroom("/analysis/lixf/sgRNA/human/LINE_screening_combine/activator.csv")
# screen_result_list <- list( L1_sup = L1_sup,
#                             L1_act = L1_act,
#                             ERV_sup = ERV_sup,
#                             ERV_act = ERV_act,
#                             SVA_sup = SVA_sup,
#                             SVA_act = SVA_act,
#                             Alu_sup = Alu_sup,
#                             Alu_act = Alu_act
#                           )

########## ENCODE RNA-seq data ############
# discard <- c("MER105","MER96","MER96B")
# HepG2df <- vroom::vroom("/analysis/lixf/RNA/human/HepG2/ENCODE/7-combine/Repeatmasker/combine_repeatmasker_DESeq2.result.tsv") %>%
#   filter(!str_detect(family,"Simple_repeat")) %>%
#   mutate(uniqueID=paste0(sample,"_",accession)) %>%
#   filter(!str_detect(subfamily,"\\?")) %>%
#   filter(!rowname %in% discard)
# K562df <- vroom::vroom("/analysis/lixf/RNA/human/K562/ENCODE/7-combine/Repeatmasker/combine_repeatmasker_DESeq2.result.tsv") %>%
#   filter(!str_detect(family,"Simple_repeat")) %>%
#   mutate(uniqueID=paste0(sample,"_",accession)) %>%
#   filter(!str_detect(subfamily,"\\?")) %>%
#   filter(!rowname %in% discard)




# human anno --------------------------------------------------------------------------------------------------------------------------


###  repeats family data
# L1=read.table('/LiuLab/reference/Human/GRCh38/TE/Repeat_Masker_hg38_sorted.L1.name',header=F)
# L1=as.vector(L1[,1])
# L1=as.character(L1)
# ERVK=read.table('/LiuLab/reference/Human/GRCh38/TE/ERVK_class.name',header=F)
# ERVK=as.vector(ERVK[,1])
# ERVK=as.character(ERVK)
# repeats=read.table('/LiuLab/reference/Human/GRCh38/TE/Repeat_Masker_hg38.NoSimpleRepeats.name',header = T)
# repeats= as.vector(repeats[,1])
# repeats=as.character(repeats)
# family <- read_tsv("/LiuLab/reference/Human/GRCh38/TE/Repeat_Masker_hg38.class",header = F)
# LTR <- family %>% filter(str_detect(V2,"LTR"))
# LTR <- as.character(LTR$X1)

# #### individual repeats data
# L1=vroom::vroom('/LiuLab/reference/Human/GRCh38/TE/TE_bed/L1.name',delim = " ",col_names =F)
# L1=as.character(L1$X1)
# ERV1=vroom::vroom('/LiuLab/reference/Human/GRCh38/TE/TE_bed/ERV1.name',delim = " ",col_names =F)
# ERV1=as.character(ERV1$X1)
# ERVK=vroom::vroom('/LiuLab/reference/Human/GRCh38/TE/TE_bed/ERVK.name',delim = " ",col_names =F)
# ERVK=as.character(ERVK$X1)
# ERVL=vroom::vroom('/LiuLab/reference/Human/GRCh38/TE/TE_bed/ERVL.name',delim = " ",col_names =F)
# ERVL=as.character(ERVL$X1)
# ERVL_MaLR=vroom::vroom('/LiuLab/reference/Human/GRCh38/TE/TE_bed/ERVL-MaLR.name',delim = " ",col_names =F)
# ERVL_MaLR=as.character(ERVL_MaLR$X1)
# SVA=vroom::vroom('/LiuLab/reference/Human/GRCh38/TE/TE_bed/SVA.name',col_names =F)
# SVA=as.character(SVA$V1)
# #### all repeats data
# repeats <- vroom::vroom("/LiuLab/reference/Human/GRCh38/TE/Repeat_Masker_hg38_sorted.bed",col_names = F)
# head(repeats)
# annotation <- vroom::vroom("/LiuLab/reference/Human/GRCh38/TE/TE_annotation.tsv",col_names = T)


#  mouse anno -------------------------------------------------------------------------------------------------------------------------
####  repeats family data
# L1=read.table('/LiuLab/reference/Mouse/GRCm38/TE/TE_class/L1.name',header=F)
# L1=as.vector(L1[,1])
# L1=as.character(L1)
# ERVK=read.table('/LiuLab/reference/Mouse/GRCm38/TE/TE_class/ERVK.name',header=F)
# ERVK=as.vector(ERVK[,1])
# ERVK=as.character(ERVK)
# repeats=read.table('/LiuLab/reference/Mouse/GRCm38/TE/Repeat_Masker_mm38.NoSimpleRepeats.name',header = T) 
# repeats= as.vector(repeats[,1]) 
# repeats=as.character(repeats)
# family <- read_tsv("/LiuLab/reference/Mouse/GRCm38/TE/TE_class/GRCm38.repeatMasker.UCSC.rmsk.name",skip = 1,col_names=F)

#### individual repeats
# L1=read.table('/LiuLab/reference/Mouse/GRCm38/TE/TE_class/L1.unique.name',header=F)
# L1=as.vector(L1[,1])
# L1=as.character(L1)
# ERV1=read.table('/LiuLab/reference/Mouse/GRCm38/TE/TE_bed/ERV1.name',header=F)
# ERV1=as.character(ERV1$V1)
# ERVL=read.table('/LiuLab/reference/Mouse/GRCm38/TE/TE_bed/ERVL.name',header=F)
# ERVL=as.character(ERVL$V1)
# ERVK=read.table('/LiuLab/reference/Mouse/GRCm38/TE/TE_bed/ERVK.name',header=F)
# ERVK=as.character(ERVK$V1)
# SVA=read.table('/LiuLab/reference/Human/GRCh38/TE/TE_bed/SVA.name',header=F)
# SVA=as.character(SVA$V1)
# ERVL_MaLR=vroom::vroom('/LiuLab/reference/Mouse/GRCm38/TE/TE_bed/ERVL-MaLR.name',delim = " ",col_names =F)
# ERVL_MaLR=as.character(ERVL_MaLR$X1)
# repeats <- vroom::vroom("/LiuLab/reference/Mouse/GRCm38/TE/GRCm38.repeatMasker.UCSC.unique.class.family.bed",col_names = F)
# head(repeats)
# annotation <- vroom::vroom("/LiuLab/reference/Mouse/GRCm38/TE/TE_annotation.tsv",col_names = T)
###### write_tsv(as.data.frame(annotation %>% filter(X7=="ERVL-MaLR") %$% anno),"/LiuLab/reference/Mouse/GRCm38/TE/TE_bed/ERVL-MaLR.name",col_names = F)
###### write_tsv(as.data.frame(annotation %>% filter(X7=="ERVK") %$% anno),"/LiuLab/reference/Mouse/GRCm38/TE/TE_bed/ERVK.name",col_names = F)
###### write_tsv(as.data.frame(annotation %>% filter(X7=="ERV1") %$% anno),"/LiuLab/reference/Mouse/GRCm38/TE/TE_bed/ERV1.name",col_names = F)
######  write_tsv(as.data.frame(annotation %>% filter(X7=="ERVL") %$% anno),"/LiuLab/reference/Mouse/GRCm38/TE/TE_bed/ERVL.name",col_names = F)


read_and_merge <- function(dir=dir,pattern="tsv",sep="\t"){
  FILE<-dir(dir,pattern=pattern,full.names = T)
  FILE <- as.data.frame(FILE)
  counts <- list()
  for (i in seq(1,nrow(FILE),by=1)) {
    counts[[i]]<- read_delim(as.character(FILE[i,]), delim = sep, col_names = TRUE,comment="#")
  }
  # combine 
  Reduce(dplyr::left_join,counts)
}

top_bottom <- function(x,n,wt,with_ties = F){
  x1 = x %>% slice_max(n=n, order_by ={{wt}},with_ties = with_ties)
  x2 = x %>% slice_min(n=n, order_by = {{wt}},with_ties = with_ties)
  x = bind_rows(x1, x2)
  return(x)
}

up_down_calculate <- function(a,p="pvalue",pcutoff=0.05,species="hg"){
  
  if(species=="hg"){  
    L1 <- a %>% dplyr::filter(str_detect(TE_ID,"L1")) %>% dplyr::filter(get(p) < pcutoff)
    #  new column indicates up and down regulated
    b <- c(rep("down",table(L1$LogFC>0)[1]),rep("up",table(L1$LogFC>0)[2]))
    #  add new column
    L1 <- L1 %>% arrange(LogFC) %>% add_column(group=b)
    L1$group <- factor(L1$group,levels=c("down","up"))
    # up and down barplot , 4 groups, "L1","L1M","L1P","FL_L1"
    # calculate number of up and down elements 
    L1_total <- L1  %>% group_by(group) %>% summarise(count=n()) %>% mutate(LINE1=rep("L1",length(group)))
    L1M <- L1 %>% dplyr::filter(str_detect(TE_ID,"L1M")) %>% group_by(group) %>% summarise(count=n()) %>% mutate(LINE1=rep("L1M",length(group)))
    L1P <- L1 %>% dplyr::filter(str_detect(TE_ID,"L1P")) %>% group_by(group) %>% summarise(count=n()) %>% mutate(LINE1=rep("L1P",length(group)))
    FLL1 <- L1 %>%  dplyr::filter(end-start > 6000) %>% group_by(group) %>% summarise(count=n()) %>% mutate(LINE1=rep("FL_L1",length(group)))
    if(length(L1_total$group)!=2){
      ifelse(L1_total$group=="up",L1_total <- L1_total %>% add_row(group="down",count=0,LINE1="L1"),L1_total <- L1_total %>% add_row(group="up",count=0,LINE1="L1"))
    }
    if(length(L1M$group)!=2){
      ifelse(L1M$group=="up",L1M <- L1M %>% add_row(group="down",count=0,LINE1="L1M"),L1M <- L1M %>% add_row(group="up",count=0,LINE1="L1M"))
    }
    if(length(L1P$group)!=2){
      ifelse(L1P$group=="up",L1P <- L1P %>% add_row(group="down",count=0,LINE1="L1P"),L1P <- L1P %>% add_row(group="up",count=0,LINE1="L1P"))
    }
    if(length(FLL1$group)!=2){
      ifelse(FLL1$group=="up",FLL1 <- FLL1 %>% add_row(group="down",count=0,LINE1="FL_L1"),FLL1 <- FLL1 %>% add_row(group="up",count=0,LINE1="FL_L1"))
    }
    # combine thme 
    c <- rbind(L1_total,L1M,L1P,FLL1)
    # convert to factor
    c$LINE1 <- factor(c$LINE1,levels=c("L1","L1M","L1P","FL_L1"))
    
  }
  
  if(species=="mm"){  
    L1 <- a %>% dplyr::filter(str_detect(TE_ID,"L1")) %>% dplyr::filter(get(p) < pcutoff)
    #  new column indicates up and down regulated
    b <- c(rep("down",table(L1$LogFC>0)[1]),rep("up",table(L1$LogFC>0)[2]))
    #  add new column
    L1 <- L1 %>% arrange(LogFC) %>% add_column(group=b)
    L1$group <- factor(L1$group,levels=c("down","up"))
    # up and down barplot , 4 groups, "L1","L1M","L1P","FL_L1"
    # calculate number of up and down elements 
    L1_total <- L1  %>% group_by(group) %>% summarise(count=n()) %>% mutate(LINE1=rep("L1",length(group)))
    L1M <- L1 %>% dplyr::filter(str_detect(TE_ID,"L1M[A-Z0-9]|HAL1|L1_Mur|L1P|L1_Rod")) %>% group_by(group) %>% summarise(count=n()) %>% mutate(LINE1=rep("Older",length(group)))
    L1Md <- L1 %>% dplyr::filter(str_detect(TE_ID,"L1Md|L1_Mus")) %>% group_by(group) %>% summarise(count=n()) %>% mutate(LINE1=rep("Younger",length(group)))
    FLL1 <- L1 %>%  dplyr::filter(end-start > 6000) %>% group_by(group) %>% summarise(count=n()) %>% mutate(LINE1=rep("FL_L1",length(group)))
    if(length(L1_total$group)!=2){
      ifelse(L1_total$group=="up",L1_total <- L1_total %>% add_row(group="down",count=0,LINE1="L1"),L1_total <- L1_total %>% add_row(group="up",count=0,LINE1="L1"))
    }
    if(length(L1M$group)!=2){
      ifelse(L1M$group=="up",L1M <- L1M %>% add_row(group="down",count=0,LINE1="Older"),L1M <- L1M %>% add_row(group="up",count=0,LINE1="Older"))
    }
    if(length(L1Md$group)!=2){
      ifelse(L1Md$group=="up",L1Md <- L1Md %>% add_row(group="down",count=0,LINE1="Younger"),L1Md <- L1Md %>% add_row(group="up",count=0,LINE1="Younger"))
    }
    if(length(FLL1$group)!=2){
      ifelse(FLL1$group=="up",FLL1 <- FLL1 %>% add_row(group="down",count=0,LINE1="FL_L1"),FLL1 <- FLL1 %>% add_row(group="up",count=0,LINE1="FL_L1"))
    }
    # combine thme 
    c <- rbind(L1_total,L1M,L1Md,FLL1)
    # convert to factor
    c$LINE1 <- factor(c$LINE1,levels=c("L1","Older","Younger","FL_L1"))
    
  }
  return(c)
}
  
  
  
  
  localization_calculate <- function(a,p="pvalue",pcutoff=0.05,species="hg"){
    if(species=="hg"){ 
      L1 <- a %>% dplyr::filter(str_detect(TE_ID,"L1")) %>% dplyr::filter(get(p) < pcutoff)
      L1_a <- L1 %>%  mutate(family=rep("L1",length(TE_ID))) %>%  group_by(family,Function) %>% summarize(count=n())
      L1M <- L1 %>% dplyr::filter(str_detect(TE_ID,"L1M")) %>%mutate(family=rep("L1M",length(TE_ID)))%>%  group_by(family,Function) %>% summarize(count=n())
      L1P <- L1 %>% dplyr::filter(str_detect(TE_ID,"L1P")) %>% mutate(family=rep("L1P",length(TE_ID)))%>%  group_by(family,Function) %>% summarize(count=n())
      FLL1 <- L1 %>%  dplyr::filter(end-start > 6000) %>% mutate(family=rep("Full Length L1",length(TE_ID)))%>%  group_by(family,Function) %>% summarize(count=n())
      c <- rbind(L1_a,L1M,L1P,FLL1)
      # convert long data to wide data
      d <- c %>% pivot_wider(names_from = Function, values_from = count)
      d <- as.data.frame(d)
      # na to 0
      d[is.na(d)] <- 0
      rownames(d) <- d$family
      d <- d[,-1]
      # calculate percentage
      return(d)
    
      }
    if(species=="mm"){
      L1 <- a %>% dplyr::filter(str_detect(TE_ID,"L1")) %>% dplyr::filter(get(p) < pcutoff)
      L1_a <- L1 %>%  mutate(family=rep("L1",length(TE_ID))) %>%  group_by(family,Function) %>% summarize(count=n())
      L1M <- L1 %>% dplyr::filter(str_detect(TE_ID,"L1M[A-Z0-9]|HAL1|L1_Mur|L1P|L1_Rod")) %>% mutate(family=rep("Older",length(TE_ID)))%>%  group_by(family,Function) %>% summarize(count=n())
      L1Md <- L1 %>% dplyr::filter(str_detect(TE_ID,"L1Md|L1_Mus")) %>% mutate(family=rep("Younger",length(TE_ID)))%>%  group_by(family,Function) %>% summarize(count=n())
      FLL1 <- L1 %>%  dplyr::filter(end-start > 6000) %>% mutate(family=rep("Full Length L1",length(TE_ID)))%>%  group_by(family,Function) %>% summarize(count=n())
      c <- rbind(L1_a,L1M,L1Md,FLL1)
      # convert long data to wide data
      d <- c %>% pivot_wider(names_from = Function, values_from = count)
      d <- as.data.frame(d)
      # na to 0
      d[is.na(d)] <- 0
      rownames(d) <- d$family
      d <- d[,-1]
      # calculate percentage
      return(d)
      
      }
    # e <- as.data.frame(t(apply(d,1,function(x) (x/sum(x)))))
    # e
  }



violin_plot <- function(a,pcutoff=0.05){
L1 <- a %>% dplyr::filter(str_detect(TE_ID,"L1")) %>% dplyr::filter(pvalue<0.05) %>% mutate(up_down=ifelse(LogFC>0,"up","down"))
sample_size = L1 %>% group_by(up_down) %>% summarize(num=n())
# Plot
L1 %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(up_down, "\n", "n=", num)) %>% 
  ggplot( aes(x=myaxis, y=abs(LogFC), fill=up_down)) +
  geom_violin(width=0.7) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2,notch = F) +
  geom_jitter(width = 0.3)+
  scale_fill_aaas() +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("") +ylab("Log2 Fold Change")+
  xlab("")+theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.title.x = element_text( size=rel(1.8)),
                 axis.title.y = element_text( size=rel(1.8)),
                 plot.title = element_text(hjust = 0.5,size=rel(2),face = "bold"),
                 axis.text= element_text(size=rel(1.5)),
                 legend.text=element_text(size=rel(1.5)),
                 legend.title = element_text(size=rel(1.5))
  )
}


plot_individual <- function(DESeqdf,class="L1",p="pvalue",pcutoff=0.05,filterExpression=200,title=NULL){
  df <- as.data.frame(DESeqdf)
  df$padj[is.na(df$padj)] <- 1
  colnames(df)[1] <- "symbol"
  
  aa <- ifelse(df$symbol %in% get(class) & df[,colnames(df)==p]< pcutoff,paste0("significant ",class),ifelse(df$symbol %in% get(class)  ,class,paste0("Non ",class)))
  
  aa[is.na(aa)] <- paste0("Non ",class)
  
  data <- df %>% add_column(group = aa)
  data$group <- factor(data$group,levels=c(paste0("Non ",class),class,paste0("significant ",class)))
  ### label data frame
  labelp <- data %>% dplyr::filter(group==paste0("significant ",class) ) %>% arrange(padj) %>% head(5)
  labelfold1<-  data %>% dplyr::filter(group==paste0("significant ",class) ) %>% top_n( 5,log2FoldChange)
  labelfold2 <-  data %>% dplyr::filter(group==paste0("significant ",class) ) %>% top_n(-5,log2FoldChange)
  labelexpression <- data %>% dplyr::filter(group==paste0("significant ",class) ) %>% top_n(5,baseMean)
  labeldf <- purrr::reduce(list(labelp,labelfold1,labelfold2,labelexpression),full_join)
  
  
  ######### unique MA plot ---------
  ggplot(data=data,mapping = aes(x=baseMean,y=log2FoldChange,color=group,size=group,alpha=group))+
    geom_point()+
    geom_hline(yintercept = 0.0, color = "red", size = 1.5,alpha=0.3)+
    scale_colour_manual(values = c("black","blue", "red"))+
    scale_alpha_manual(values = c(0.15,0.2,1))+
    scale_size_manual(values = c(0.1,0.9,1.5))+
    ggrepel::geom_text_repel(data= labeldf,label=str_split_fixed(labeldf$symbol,pattern = "_chr",2)[,1],size=7, force = 2,segment.size=0.2,segment.alpha=0.95,segment.color="black",color="black")+
    scale_x_continuous(labels = comma,trans = "log10")+
    #ylim(-2,2)+
    xlab("RNA-seq reads (normalized counts)")+
    ylab("Gene expression change (log2[Treatment/control])") +
    labs(subtitle = "RepeatMasker Unique")+
    ggtitle(paste0(title," , ",p," < ",pcutoff,", ",class))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text( size=rel(1.8)),
          axis.title.y = element_text( size=rel(1.8)),
          plot.title = element_text(hjust = 0.5,size=rel(2),face = "bold"),
          plot.subtitle = element_text(hjust = 1,size=rel(1.2)),
          axis.text.x = element_text(angle=0,vjust=0),
          axis.text= element_text(size=rel(1.5)),
          legend.text=element_text(size=rel(1.5)),
          legend.title = element_text(size=rel(1.5)),
          legend.position = "bottom",
          #plot.margin = margin(0.4,1,0.3,0.3, "cm"),
          strip.text.x = element_text(size = rel(1.8)),
          strip.placement = "outside",strip.background = element_blank()
    )+
    guides(color=guide_legend(title="",override.aes = list(size=5)),alpha="none",size="none")
}


get_top_label <- function(data,group,n){
  labelp <- data %>% dplyr::filter(group=={{group}} ) %>% arrange(padj) %>% head(n)
  labelfold1<-  data %>% dplyr::filter(group=={{group}}) %>% top_n( n,log2FoldChange)
  labelfold2 <-  data %>% dplyr::filter(group=={{group}}) %>% top_n(-n,log2FoldChange)
  labelexpression <- data %>% dplyr::filter(group=={{group}} ) %>% top_n(n,baseMean)
  labeldf <- purrr::reduce(list(labelp,labelfold1,labelfold2,labelexpression),full_join)
  return(labeldf)
}
### family <- read_tsv("/LiuLab/reference/Human/GRCh38/TE/Repeat_Masker_hg38.class",col_names=F)
RepeatMAplot <- function(DESeqdf,p="pvalue",pcutoff=0.05,class="L1",filterExpression=0,title="",label_top_n=5){
  df <- as.data.frame(DESeqdf)
  df$padj[is.na(df$padj)] <- 1
  colnames(df)[1] <- "symbol"
  if (class=="family") {
    ERV1 <- family %>% filter(str_detect(V3,"ERV1"))
    ERV1 <- as.character(ERV1$V1)
    ERVK <- family %>% filter(str_detect(V3,"ERVK"))
    ERVK <- as.character(ERVK$V1)
    ERVL <- family %>% filter(str_detect(V3,"ERVL$"))
    ERVL <- as.character(ERVL$V1)
    ERVL_MaLR <- family %>% filter(str_detect(V3,"ERVL-MaLR"))
    ERVL_MaLR <- as.character(ERVL_MaLR$V1)
    LINE1 <- family %>% filter(str_detect(V3,"L1"))
    LINE1 <- as.character(LINE1$V1)
    
    aa <- case_when(df$symbol %in% ERV1 & df[,colnames(df)==p] < pcutoff ~ "sig ERV1",
                    df$symbol %in% ERVL & df[,colnames(df)==p] < pcutoff ~ "sig ERVL",
                    df$symbol %in% ERVK & df[,colnames(df)==p] < pcutoff ~ "sig ERVK",
                    df$symbol %in% ERVL_MaLR & df[,colnames(df)==p] < pcutoff ~ "sig ERVL_MaLR",
                    df$symbol %in% LINE1 & df[,colnames(df)==p] < pcutoff ~ "sig LINE1",
                    TRUE ~ "others"
                    )
    aa[is.na(aa)] <- paste0("others")
    data <- df %>% add_column(group = aa)
    
    data$group <- factor(data$group,levels=c("others" ,"sig LINE1", "sig ERV1","sig ERVK","sig ERVL","sig ERVL_MaLR" ))
    
    tmp <- unique(data$group)
    ### label data frame
    combine_list <- lapply(tmp[tmp!="others"],function(x) get_top_label(data=data,group=x,label_top_n))
    ### test list empty
    if(length(combine_list)==0){
      # select empty df
      labeldf <- data %>% filter(symbol=="random")
    } else{
      labeldf <- purrr::reduce(combine_list,full_join)
    }
    #####
    plot <- ggplot(data=data,mapping = aes(x=baseMean,y=log2FoldChange,color=group,size=group,alpha=group))+
      geom_point()+
      geom_hline(yintercept = 0.0, color = "red", size = 1.5,alpha=0.3)+
      scale_colour_manual(values = c("grey",'#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'))+
      scale_alpha_manual(values = c(0.5,1,1,1,1,1))+
      scale_size_manual(values = c(0.5,2.5,2.5,2.5,2.5,2.5))+
      ggrepel::geom_text_repel(data=labeldf,label=labeldf$symbol,size=7, force = 2,segment.size=0.2,segment.alpha=0.95,segment.color="black",color="black")+
      scale_x_continuous(labels = comma,trans = "log10")+
      #ylim(-1.5,1.5)+
      xlab("RNA-seq reads (normalized counts)")+
      ylab("Gene expression change (log2[Treatment/control])") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.x = element_text( size=rel(1.2)),
            axis.title.y = element_text( size=rel(1.2)),
            plot.title = element_text(hjust = 0.5,size=rel(1.5),face = "bold"),
            plot.subtitle = element_text(hjust = 1,size=rel(1.2)),
            axis.text.x = element_text(angle=0,vjust=0),
            axis.text= element_text(size=rel(1.2)),
            legend.text=element_text(size=rel(1.2)),
            legend.title = element_text(size=rel(1.2)),
            legend.position = "bottom",
            plot.margin = margin(0.4,1,0.3,0.3, "cm"),
            strip.text.x = element_text(size = rel(1.5)),
            strip.placement = "outside",strip.background = element_blank()
      )+
      guides(color=guide_legend(title="",override.aes = list(size=3)),alpha="none",size="none")+
      ggtitle(paste0(title," , ",p," < ",pcutoff,", ",class))
    plot
  } else{
    aa <- ifelse(df$symbol %in% get(class) & df[,colnames(df)==p]< pcutoff,paste0("significant ",class),ifelse(df$symbol %in% get(class)  ,class,paste0("others")))
    aa[is.na(aa)] <- paste0("others")
    data <- df %>% add_column(group = aa)
    data$group <- factor(data$group,levels=c(paste0("others"),class,paste0("significant ",class)))
    ### label data frame
    labeldf <- get_top_label(data=data,group=paste0("significant ",class),label_top_n)
    
    plot <- ggplot(data=data,mapping = aes(x=baseMean,y=log2FoldChange,color=group,size=group,alpha=group))+
      geom_point()+
      geom_hline(yintercept = 0.0, color = "red", size = 1.5,alpha=0.3)+
      scale_colour_manual(values = c("grey","blue", "red"))+
      scale_alpha_manual(values = c(0.5,0.7,1))+
      scale_size_manual(values = c(0.5,1.5,2.5))+
      ggrepel::geom_text_repel(data=labeldf,label=labeldf$symbol,size=7, force = 2,segment.size=0.2,segment.alpha=0.95,segment.color="black",color="black")+
      scale_x_continuous(labels = comma,trans = "log10")+
      #ylim(-1.5,1.5)+
      xlab("RNA-seq reads (normalized counts)")+
      ylab("Gene expression change (log2[Treatment/control])") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.x = element_text( size=rel(1.8)),
            axis.title.y = element_text( size=rel(1.8)),
            plot.title = element_text(hjust = 0.5,size=rel(2),face = "bold"),
            plot.subtitle = element_text(hjust = 1,size=rel(1.2)),
            axis.text.x = element_text(angle=0,vjust=0),
            axis.text= element_text(size=rel(1.5)),
            legend.text=element_text(size=rel(1.5)),
            legend.title = element_text(size=rel(1.5)),
            legend.position = "bottom",
            plot.margin = margin(0.4,1,0.3,0.3, "cm"),
            strip.text.x = element_text(size = rel(1.8)),
            strip.placement = "outside",strip.background = element_blank()
      )+
      guides(color=guide_legend(title="",override.aes = list(size=5)),alpha="none",size="none")+
      ggtitle(paste0(title," , ",p," < ",pcutoff,", ",class))
  }
  #p2 <- ggMarginal(p,type = "density",margins = "y",groupColour  = T,groupFill = T)
  return(plot)
}

individual_repeats_boxplot <- function(DESeqdf,annotation=annotation,sig=FALSE,p="pvalue",pcutoff=0.05,class="L1",title=""){
  df <- as.data.frame(DESeqdf)
  df$padj[is.na(df$padj)] <- 1
  colnames(df)[1] <- "symbol"
  
  #aa <- ifelse(df$symbol %in% get(class) & df[,colnames(df)==p]< pcutoff,paste0("significant ",class),ifelse(df$symbol %in% get(class)  ,class,paste0("Non ",class)))
  #aa <- ifelse(df$symbol %in% L1 & df$pvalue<0.05,"significant L1",ifelse(df$symbol %in% L1  ,"LINE-1","Non LINE-1")) 
  if(sig==TRUE){
    data <- df %>%  dplyr::filter( .data[[p]] <= pcutoff)#%>% add_column(group = aa)
  }
  else  {
    data <- df 
  }
  repeatclass <- left_join(data %>% mutate(anno=symbol),
                     annotation %>% dplyr::select(anno,X5,X6,X7),
                     by="anno")
  #table(repeatclass$X7)
  # tmp <- repeatclass %>% mutate(repeatclass=ifelse(group==paste0("significant ",class), paste0("significant ",class),X7))
  tmp <- repeatclass %>% mutate(repeatclass=X7)
  #### filter low number repeats
  names(sort(table(tmp$repeatclass),decreasing = T))
  list1 <- c("Alu","Simple_repeat","L1","MIR","L2","ERV1","Low_complexity","ERVL-MaLR","hAT-Charlie","ERVL","TcMar-Tigger","CR1","ERVK","significant L1","SVA","Satellite")
  list2 <- c("Simple_repeat","Low_complexity","Alu","B2","ERVK","L1","B4","ERVL-MaLR","Satellite","MIR","L2","hAT-Charlie","ID","ERVL","ERV1","CR1","tRNA","TcMar-Tigger","rRNA","significant L1")
  list <- c(list1,list2)
  tmp1 <- tmp %>% filter(repeatclass %in% list)
  table(tmp1$repeatclass)
  ### create name with number
  numberdf <- tmp1 %>% group_by(repeatclass) %>% summarize(count=n(),median=median(log2FoldChange)) %>% mutate(group2=paste0(repeatclass,"\n","(",count,")"))
  # get order by median log2FC
  levels <- arrange(numberdf,median) %$% group2
  # data frame used for plotting, combine number data frmae
  plotdf <- left_join(tmp1, numberdf %>% dplyr::select(repeatclass,group2))
  # order of plotting
  plotdf$group2  <- factor(plotdf$group2,levels=levels)
  ####  used for multiple comparasion
  # group2 <- as.character(unique(plotdf$group2))
  # my_comparation <- list(group2[c(1,2)],group2[c(1,3)],group2[c(2,3)])
  
  p2 <- ggplot(plotdf,aes(x=group2,y=log2FoldChange,fill="#group2"))+
    geom_boxplot(width=0.6)+
    geom_hline(yintercept = 0,color="red",alpha=0.5,size=1.2)+
    coord_flip()+
    ###  multiple comparasion
    #  stat_compare_means(comparisons = my_comparation,size=6,label = "p.signif")+
    #  stat_compare_means(label.y = 5)+
    # scale_fill_manual(values=c("#00A087FF","#4DBBD5FF","#E64B35FF"))+
    scale_fill_manual(values=c("#4DBBD5FF"))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text( size=rel(1.8)),
          axis.title.y = element_text( size=rel(1.8)),
          plot.title = element_text(hjust = 0.5,size=rel(2),face = "bold"),
          plot.subtitle = element_text(hjust = 1,size=rel(1.2)),
          axis.text.x = element_text(angle=0,vjust=0),
          axis.text= element_text(size=rel(1.5)),
          legend.text=element_text(size=rel(1.5)),
          legend.title = element_text(size=rel(1.5)),
          plot.margin = margin(0.3,1,0.3,0.3, "cm"),
          strip.text.x = element_text(size = rel(1.8)),
          legend.position = "none",
          strip.placement = "outside",strip.background = element_blank()
    )+
    guides(fill=guide_legend(title=""))+
    xlab("")
  if(sig==TRUE){
   p2 <- p2+ggtitle(paste0(title," , ",p," < ",pcutoff))
  }
  else  {
   p2 <- p2+ggtitle(title)
  }
  p2
  }



countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToEffCounts <- function(counts, len, effLen)
{
  counts * (len / effLen)
}

dft <- function(c) {
  mtotal <- c[,-1]
  mtotal <- t(as.matrix(mtotal))
  colnames(mtotal) <- unlist(c[,1])
  return(mtotal)
  
  # x <- t(c)
  # colnames(x) <- x[1,]
  # x <- x[-1,]
  # y <- apply(x,2,as.numeric)
  # rownames(y) <- rownames(x)
  # y <- as.data.frame(x) %>% rownames_to_column()
  # y
}

convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "dec2021.archive.ensembl.org")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  
  # Print the first 6 genes found to the screen
  
  return(genesV2)
}

#library(dplyr)

convert_mouse_to_human <- function(gene_list){
  mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  
  return (output)
}



plot_cNMF <- function(dir,gene_spectra,usage,count_input,topgene=100,metadata,cutoff_value=0.2){
  setwd(paste0(dir))
  
  col57 <- c('#9af764', '#3e82fc', '#fe0002', '#f4d054', '#ed0dd9',
             '#13eac9', '#e4cbff', '#b1d27b', '#ad8150', '#601ef9',
             '#ff9408', '#75bbfd', '#fdb0c0', '#a50055', '#4da409',
             '#c04e01', '#d2bd0a', '#ada587', '#0504aa', '#650021',
             '#d0fefe', '#a8ff04', '#fe46a5', '#bc13fe', '#fdff52',
             '#f2ab15', '#fd4659', '#ff724c', '#cba560', '#cbf85f',
             '#78d1b6', '#9d0216', '#874c62', '#8b88f8', '#05472a',
             '#b17261', '#a4be5c', '#742802', '#3e82fc', '#eedc5b',
             '#a8a495', '#fffe71', '#c1c6fc', '#b17261', '#ff5b00',
             '#f10c45', '#3e82fc', '#de9dac', '#f10c45', '#056eee',
             '#e6daa6', '#eedc5b', '#c87606', '#9dbcd4', '#56ae57',
             '#49759c', '#d8dcd6')
  # dir="/analysis2/lixf/lulab_data/total/10-NMF/up_intron_inter_nosimple_HCC/"
  # metadata <- fread("/analysis2/lixf/lulab_data/metadata.txt")
  # spectra <- fread(paste0(dir,"/up_intron_inter_nosimple_HCC.gene_spectra_score.k_2.dt_2_00.txt"), header = TRUE)
  # usage <- fread(paste0(dir,"/up_intron_inter_nosimple_HCC.usages.k_2.dt_2_00.consensus.txt"),header=TRUE)
  # total <- fread("/analysis2/lixf/lulab_data/total/10-NMF/up_intron_inter_nosimple_HCC.total.cNMF.input.txt")
  
  #total[1:5,1:5]
  mtotal <- count_input[,-1]
  mtotal <- t(as.matrix(mtotal))
  colnames(mtotal) <- as.vector(t(count_input[,1]))
  #mtotal[1:5,1:5]
  
  spectra <- t(gene_spectra[,-1])
  
  topgenes <- list()
  ###  number of cluster
  cluster_num <- ncol(spectra)
  ## get the top 100 genes in each cluster
  for(i in 1:cluster_num){
    topgenes[[i]] <- names(sort(spectra[,i], decreasing = TRUE))[1:topgene]
  }
  names(topgenes) <- c(paste0("cluster",c(1:cluster_num)))
  
  ## define the cluster for all samples
  usage <- as.data.frame(usage)
  rownames(usage) <- as.vector(t(usage[,1]))
  usage <- usage[,-1]
  colnames(usage) <- paste0("cluster",1:cluster_num)
  usage_norm <- usage/rowSums(usage)
  #head(usage_norm)
  
  usage_class <- data.frame(usage_norm,
                            maxUsage = apply(usage_norm, 1, function(x) max(x)),
                            max2 = apply(usage_norm, 1, function(x) sort(x, decreasing = TRUE)[2]),
                            cluster = apply(usage_norm, 1, function(x) colnames(usage_norm)[which.max(x)]))
  
  usage_class <- data.frame(usage_class, diff = usage_class$maxUsage-usage_class$max2)
  
  #head(usage_class)
  #table(usage_class$maxUsage > 0.3)
  #table(usage_class$diff > 1/6)
  #table(usage_class$cluster)
  
  usage_cluster <- list()
  for(i in 1:cluster_num){
    usage_cluster[[i]] <- rownames(usage_class)[usage_class$cluster%in%c(paste0("cluster",i))]
  }
  #### plot the top genes - samples heatmap
  #mtotal[unlist(topgenes),unlist(usage_cluster)]
  #heatmap_mat <- RNAheatmap_Plotmat(t(mtotal[unlist(topgenes),unlist(usage_cluster)]), featuregene = topgenes, samplesCluster = usage_cluster, samplekeep = 1, genekeep = 0)
  # mtotal
  #
  # genelist <- NULL
  # samplelist <- NULL
  # featuregene = topgenes; samplesCluster = usage_cluster
  heatmap_mat <- mtotal[unlist(topgenes),unlist(usage_cluster)]
  dim(heatmap_mat)
  heatmap_mat[1:4,1:4]
  cluster_num
  unlist(lapply(c(1:cluster_num),function(x) rep(paste0("clust",x),length(topgenes[[x]]))))
  
  
  # rowanno <-data.frame(FeatureTE = c(rep("clust1", length(topgenes$cluster1)), rep("clust2", length(topgenes$cluster2)),
  #                                    rep("clust3", length(topgenes$cluster3)), rep("clust4", length(topgenes$cluster4)),
  #                                    rep("clust5", length(topgenes$cluster5)), rep("clust6", length(topgenes$cluster6))))
  rowanno <- data.frame(FeatureTE = unlist(lapply(c(1:cluster_num),function(x) rep(paste0("clust",x),length(topgenes[[x]])))))
  
  ###rownames(rowanno) <- make.names(unlist(topgenes), unique = TRUE)
  rownames(rowanno) <- unlist(topgenes)
  
  
  genes <- col57[30:(30+cluster_num-1)]
  names(genes) <- names(table(rowanno$FeatureTE))
  
  ## annotation information for each sample
  metadata <- metadata %>% as.data.frame()
  rownames(metadata) <- metadata$sample_id
  colanno <- data.frame(cluster = usage_class[colnames(heatmap_mat),"cluster"], metadata[colnames(heatmap_mat), c("source","labels","stage")])
  colnames(colanno)[c(2,3)] <- c("hospitals","cancers")
  rownames(colanno) <- colnames(heatmap_mat)
  #head(colanno)
  colanno[is.na(colanno)] <- "NA"
  
  #table(colanno$cluster) #cluster1 cluster2 cluster3 cluster4 cluster5 cluster6
  #table(colanno$hospitals) # I   IA   II  IIA  IIB  IIC  III IIIA IIIB IIIC   IV  IVA  IVB   NA
  #length(table(colanno$cancers)) # MSI-H MSI-L   MSS    NA
  #length(table(colanno$stage)) # CIN     GS HM-SNV    MSI     NA
  
  cluster <- col57[1:(1+cluster_num-1)]
  names(cluster) <- names(table(colanno$cluster))
  #
  #metadata
  hospitals <- c(col57[6:12])
  names(hospitals) <- names(table(colanno$hospitals))
  
  cancers <- c(col57[21:26])
  names(cancers) <- names(table(colanno$cancers))
  
  stage <- c(col57[24:35])
  names(stage) <- names(table(colanno$stage))
  
  require(pheatmap)
  #tmp <- heatmap_mat[1:4,1:4]
  #scale(t(tmp))
  ## scale the expression matrix by gene
  tmp_mat <- scale(t(heatmap_mat))
  tmp_mat[1:4,1:4]
  #tmp_mat[(is.na(tmp_mat))] <- 0
  
  quantile(tmp_mat)
  heatmap_mat[1:4,1:4]
  #tmp_mat <- t(tmp_mat)
  
  tmp_mat[tmp_mat > cutoff_value] <- cutoff_value
  tmp_mat[tmp_mat < -cutoff_value] <- -cutoff_value
  
  
  p <- pheatmap(t(tmp_mat),
                #scale = "row",
                #color = c(colorRampPalette(c("green", "black"))(n = 100),colorRampPalette(c("black", "red"))(n = 100)),
                color = c(colorRampPalette(c("#00008c", "white"))(n = 100),colorRampPalette(c("white", "#7E160E"))(n = 100)),
                #color = colorRampPalette(c("white", "#7E160E"))(n = 100),
                #color = colorRampPalette(c("darkgreen", "darkred"))(n = 100),
                #color = c(colorRampPalette(c("blue", "white"))(n = 100),colorRampPalette(c("white", "red"))(n = 100)),
                
                border_color = "white",
                #legend = FALSE,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                labels_col = FALSE,
                annotation_col = colanno,
                annotation_row = rowanno,
                show_rownames = FALSE,
                annotation_colors = list(genes = genes,cluster=cluster,stage = stage, hospitals = hospitals,cancers = cancers, cluster = cluster)
  )
  
  list <- list(heatmap=p,
               topgenes=topgenes,
               usage_cluster=usage_cluster,
               heatmap_mat=tmp_mat)
  return(list)
}

gene_to_mat <- function(knownTPMdf,genelist){
  genelist_known <- knownTPMdf %>% filter(GeneSymbol %in% genelist)
  genelist_known$GeneSymbol <- factor( genelist_known$GeneSymbol,levels=genelist)
  genelist_known <- genelist_known %>% arrange(GeneSymbol)
  
  mat <- genelist_known[,-1]
  mat <- apply(mat,2,function(x) log2(x+1) )
  rownames(mat) <- genelist_known$GeneSymbol
  return(mat)
}



# stemness calculate ------------------------------------------------------------------------------------------------------------------
### TPMdf example, row as gene, column as sample, return each sample stemness score
## stem <-  calculate_stemness(TPMdf,fnOut = "NCCIT_stemness.tsv")
# TPMdf[1:5,1:5]
# A tibble: 5 x 5
# GeneSymbol  KDM2B_KO_CxxC_1_combine.sorted.bam KDM2B_KO_CxxC_2_combine.sorted.bam KDM2B_KO_CxxC_3_combine.sorted.bam KDM2B_KO_exon1_2_combine.sorted.bam
# <chr>                                    <dbl>                              <dbl>                              <dbl>                               <dbl>
#   1 DDX11L1                                    0                                0.137                             0.0881                                0   
# 2 WASH7P                                    14.1                             17.4                              17.0                                   8.31
# 3 MIR6859-1                                 33.4                             10.5                              42.7                                  30.4 
# 4 MIR1302-2HG                                0                                0                                 0                                     0   
# 5 MIR1302-2                                  0                                0                                 0                                     0   
calculate_stemness <- function( fnTPM=TPMdf, fnOut = "mRNAsi_stemness.tsv" )
{
  
  ## Load the signature
  #w <- read.delim( fnSig, col_names=FALSE, row.names=1 ) %>% as.matrix() %>% drop()
  w <- vroom::vroom( "/analysis/lixf/stemness/pcbc-stemsig.tsv", col_names=FALSE )  %>% 
    column_to_rownames("X1") %>% 
    as.matrix() %>% 
    drop()
  
  ## Reduces HUGO|POSITION gene IDs to just HUGO
  # f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )
  
  #s <- synGet( "syn4976369", downloadLocation = "/data/pancan" )
  # X <- vroom::vroom( "./pancan.tsv" ) %>% #Read the raw values
  # 	filter( !grepl( "\\?", gene_id ) )  %>%      ## Drop genes with no mapping to HUGO
  # 	mutate( gene_id = f( gene_id ) ) %>%        ## Clip gene ids to HUGO
  # 	filter( gene_id %in% names(w) )         ## Reduce to the signature's gene set
  # 
  X <- fnTPM %>% 
    filter(GeneSymbol %in% names(w) )
  
  ## SLC35E2 has multiple entries with the same HUGO id
  ## Keep the first entry only
  # j <- grep( "SLC35E2", X[,1] )
  # if( length(j) > 1 )
  # 	X <- X[-j[-1],]
  
  ## Convert to a matrix
  rownames(X) <- NULL
  X <- X %>% tibble::column_to_rownames( "GeneSymbol" ) %>% as.matrix()
  
  ## Reduce the signature to the common set of genes
  stopifnot( all( rownames(X) %in% names(w) ) )
  w <- w[ rownames(X) ]
  
  ####### Score via Spearman correlation
  s <- apply( X, 2, function(z) {cor( z, w, method = "sp", use = "complete.obs" )} )
  
  ## Scale the scores to be between 0 and 1
  s <- s - min(s)
  s <- s / max(s)
  
  write.table(cbind(s), file = fnOut, sep = "\t", quote = FALSE, col.names = FALSE)
  s
}

#### screnn_list like /analysis/lixf/sgRNA/human/ERV_screening_combine/suppressor.csv
require(patchwork)
require(scales)
get_anno_gene <- function(gene,size=5){
  # ERV_sup <- vroom::vroom("/analysis/lixf/sgRNA/human/ERV_screening_combine/suppressor.csv")
  # ERV_act <- vroom::vroom("/analysis/lixf/sgRNA/human/ERV_screening_combine/activator.csv")
  # L1_sup <- vroom::vroom("/analysis/lixf/sgRNA/human/LINE_screening_combine/suppressor.csv")
  # L1_act <- vroom::vroom("/analysis/lixf/sgRNA/human/LINE_screening_combine/activator.csv")
  # get target gene
  tmp <- L1_sup %>% filter(id == gene)
  rankplot <- function(rankset){
    tmp <- rankset %>% filter(id == gene) 
    tmp_melt <- reshape2::melt(tmp,id.vars = "id",measure.vars = colnames(tmp)[grep("rank",colnames(tmp))] )
    tmp_melt$variable <- str_remove(tmp_melt$variable,"\\.rank")
    tmp_melt$variable <- str_replace(tmp_melt$variable,"_plus"," +")
    tmp_melt$variable <- str_replace(tmp_melt$variable,"_minus"," -")
    tmp_melt$variable <- str_replace(tmp_melt$variable,"_1"," 1st")
    tmp_melt$variable <- str_replace(tmp_melt$variable,"_2"," 2nd")
    ### rank plot 
    ggplot(tmp_melt,aes(x=variable,y=value)) +
      geom_bar(stat="identity")+
      geom_text(aes(label=value), 
                position=position_dodge(width=0.9),vjust=-0.25,hjust=0.5,size=rel(4))+
      scale_y_continuous(trans = "log10",limits = c(1,20000))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.x = element_text( size=rel(1.8)),
            axis.title.y = element_text( size=rel(1.8)),
            plot.title = element_text(hjust = 0.5,size=rel(2),face = "bold"),
            plot.subtitle = element_text(hjust = 1,size=rel(1.2)),
            axis.text.x = element_text(angle=30,vjust=0.5),
            axis.text= element_text(size=rel(1.5)),
            legend.text=element_text(size=rel(1.5)),
            legend.title = element_text(size=rel(1.5)),
            plot.margin = margin(0.3,1,0.3,0.3, "cm"),
            strip.text.x = element_text(size = rel(1.8)),
            strip.placement = "outside",strip.background = element_blank(),
            legend.position = "bottom"
      )+
      guides(fill=guide_legend(title=""))+xlab("")+ylab("rank")+ggtitle(deparse(substitute(rankset)))
    
    }
  pERV_sup <- rankplot(ERV_sup)
  pERV_act <- rankplot(ERV_act)
  pL1_sup <- rankplot(L1_sup)
  pL1_act <- rankplot(L1_act)
 #%>% select(  id,`Subcellular Location (Protein Atlas)`,`Gene Summary`,`Biological Process (GO)`,contains("rank"),contains("desc") )

  ## GO BP
  
  GO <- str_split(tmp$`Biological Process (GO)`,";")
  GO <- (sapply(GO,function(x) str_wrap(x,width=70)))
  GO <- paste(GO, sep="", collapse="\n")
  
  ## Go MF
  
  GOMF <- str_split(tmp$`Molecular Function (GO)`,";")
  GOMF <- (sapply(GOMF,function(x) str_wrap(x,width=70)))
  GOMF <- paste(GOMF, sep="", collapse="\n")
  
  ## GO CC
  
  GOCC <- str_split(tmp$`Cellular Component (GO)`,";")
  GOCC <- (sapply(GOCC,function(x) str_wrap(x,width=70)))
  GOCC <- paste(GOCC, sep="", collapse="\n")
  
  ### localization
  
  localization <- tmp$`Subcellular Location (Protein Atlas)`
  
  ### rank plot 

  ### GO BP
  p2 <- ggplot()+
    annotate("text",
              x = 5,
              y = 5,
              size = size,
              label = GO ) +
  #xlim(c(0,10))+ylim(c(0,10))+
    ggtitle("GO Biological Process")+
    theme_void()+
    theme(plot.title = element_text(size = 15,hjust = 0,face = "bold"))
  
  ### GO MF 
  p3 <- ggplot()+
    annotate("text",
             x = 5,
             y = 5,
             size = size,
             label = GOMF ) +
    #xlim(c(0,10))+ylim(c(0,10))+
    ggtitle("GO Molecular Function")+
    theme_void()+
    theme(plot.title = element_text(size = 15,hjust = 0,face = "bold"))
  
  ## GO CC
  p4 <- ggplot()+
    annotate("text",
             x = 5,
             y = 5,
             size = size,
             label = GOCC ) +
  #xlim(c(0,10))+ylim(c(0,10))+
    ggtitle("GO Cellular Component")+
    theme_void()+
    theme(plot.title = element_text(size = 15,hjust = 0,face = "bold"))
  
  ### Subcellular Localization
  p5  <- ggplot()+
    annotate("text",
             x = 5,
             y = 7,
             size = 4,
             label = localization ) +
    #xlim(c(0,10))+ylim(c(0,10))+
    ggtitle("Subcellular Location")+
    theme_void()+
    theme(plot.title = element_text(size = 15,hjust = 0,face = "bold"))
  
  patch <- ((pERV_sup+pERV_act)/(pL1_sup+pL1_act) ) |(p2/p3/p4/p5+plot_layout(heights = c(3,3,3,1)) )
  patch <- patch+
    plot_layout(width = c(1,2))+
    plot_annotation(title=gene,
                  subtitle = tmp$Description,
                  theme = theme(plot.title = element_text(size = 25,hjust=0.5),
                                plot.subtitle = element_text(size = 15,hjust=0.5))
                  )
  patch
}





### anno bed was the output of 6-LINE1M.sh, use the XXXXXX.host.annotated.bed file.

calculate_individual_updown <- function(annobed,title="",filter_copy_number=5,p="pvalue",pcutoff=0.05,species = "hg"){
  list <- list()
  if (species == "hg") {
    annotation <- vroom::vroom("/LiuLab/reference/Human/GRCh38/TE/TE_annotation.tsv",col_names = T)  
  } 
  if (species == "mm" ) {
    annotation <- vroom::vroom("/LiuLab/reference/Mouse/GRCm38/TE/TE_annotation.tsv",col_names = T)
  }
  
  family <- annotation %>% select(X5,X6,X7) %>% distinct()
  colnames(family) <- c("subfamily","class","family")
  
  annodf <- read_tsv(annobed) %>% 
    filter(UQ(rlang::sym(p))   < pcutoff) %>% 
    filter(!str_detect(TE_ID,"\\(")) %>% 
    mutate(subfamily = str_split_fixed(TE_ID,"_chr",2)[,1] )
  anno_family <- left_join(annodf,family,by="subfamily")
  
  
  tmp <- anno_family %>% group_by(family) %>% summarize(count = n()) %>% filter(count > filter_copy_number )%$%family
  # updown ------------------------------------------------------------------------------------------------------------------------------
  updown <- anno_family %>% 
    mutate(updown = ifelse(LogFC > 0, "up","down")) %>% 
    group_by(family,updown) %>% 
    summarize(count = n()) %>%
    filter(family %in% tmp ) %>% 
    group_by(family) %>% 
    mutate(percent = count/sum(count),
           label_y = cumsum(percent) - 0.5 * percent,
           updown = factor(updown, levels = c("down","up")))
  ### control levels 
  levels_up <- updown %>% filter(updown == "up") %>% arrange(percent)%$%family
  levels_others <- updown$family[! (updown$family %in% levels_up)]
  levels <- c(levels_others,levels_up)
  updown$family <- with(updown,factor(family,levels = levels))
  
  anno_df <- updown %>% group_by(family) %>% summarize(count=sum(count)) %>% mutate(family_anno = paste0(family," (",count,")")) %>% select(family,family_anno)
  updown <- left_join(updown,anno_df,by="family")
  
  levels_up <- updown %>% filter(updown == "up") %>% arrange(percent)%$%family_anno
  levels_others <- updown$family_anno[! (updown$family_anno %in% levels_up)]
  levels <- c(levels_others,levels_up)
  updown$family_anno <- factor(updown$family_anno,levels = levels)
  
  p <- updown %>% 	
    ggplot(aes(x=family_anno,y=percent,fill=updown))+
    geom_bar(stat="identity",position="stack")+
    scale_y_continuous(labels = scales::percent)+
    ggsci::scale_fill_aaas()+
    #geom_text(aes(y = label_y, label = count), hjust = 0.5, colour = "white")+
    coord_flip()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text( size=rel(1.8)),
          axis.title.y = element_text( size=rel(1.8)),
          plot.title = element_text(hjust = 0.5,size=rel(2),face = "bold"),
          plot.subtitle = element_text(hjust = 1,size=rel(1.2)),
          axis.text.x = element_text(angle=0,vjust=0),
          axis.text= element_text(size=rel(1.2)),
          legend.text=element_text(size=rel(1.5)),
          legend.title = element_text(size=rel(1.5)),
          plot.margin = margin(0.3,1,0.3,0.3, "cm"),
          strip.text.x = element_text(size = rel(1.8)),
          strip.placement = "outside",strip.background = element_blank(),
          legend.position = "bottom"
    )+
    guides(fill=guide_legend(title=""))+
    labs(title = paste0(title),
         subtitle = "individual repeats expression change",
         caption = paste0("copy number >",filter_copy_number," ", p," < ", pcutoff, " DESeq2")
    )+
    xlab("")+ylab("")
  list <- list(updown,p)
  return(list)
}

calculate_localization <- function(annobed,class="L1",title=" ",p="pvalue",pcutoff=0.05,species = "hg"){
  list <- list()
  if (species == "hg") {
    annotation <- vroom::vroom("/LiuLab/reference/Human/GRCh38/TE/TE_annotation.tsv",col_names = T)  
  } 
  if (species == "mm" ) {
    annotation <- vroom::vroom("/LiuLab/reference/Mouse/GRCm38/TE/TE_annotation.tsv",col_names = T)
  }
  
  family <- annotation %>% select(X5,X6,X7) %>% distinct()
  colnames(family) <- c("subfamily","class","family")
  
  annodf <- read_tsv(annobed) %>% 
    filter(UQ(rlang::sym(p))   < pcutoff) %>% 
    filter(!str_detect(TE_ID,"\\(")) %>% 
    mutate(subfamily = str_split_fixed(TE_ID,"_chr",2)[,1] )
  anno_family <- left_join(annodf,family,by="subfamily")
  
  functions <- c("others","splicing","downstream","upstream" ,  "UTR3"    ,   "UTR5"  ,  "intergenic" ,"exonic" ,  "intronic"  )
  colors <- c("#808180FF","#A20056FF","#5F559BFF","#3B4992FF", "#EE0000FF" ,"#008B45FF" ,"#631879FF", "#008280FF" ,"#BB0021FF")
  
  df <- anno_family[str_detect(anno_family$TE_ID,class),]
  p <- df %>% 
    mutate(Function = str_remove(Function,"ncRNA_"),
           Function = ifelse(Function %in% functions,Function,"others")) %>% 
    group_by(Function) %>% 
    summarize(count=n()) %>%
    arrange(count) %>% 
    mutate(Function = fct_inorder(Function)) %>% 
    ggplot(aes(x="", y=count, fill=Function)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = colors,
                      breaks = functions,
                      guide = guide_legend(reverse = TRUE) )+
    ggtitle(paste0(title,class))+
    theme_void()+
    guides(fill=guide_legend(title=""))
  p
  return(p)
}


# get TCGA gene-repeat correaltion ----------------------------------------


# gene_expression <- fread("/analysis2/lixf/TCGA/corrplot/TCGA_gene_expression.csv")
# 
# repeat_expression <- fread("/analysis2/lixf/TCGA/corrplot/TCGA_repeat_expression.csv")

get_gene_repeat_corr <- function(gene_expression,repeat_expression,gene,repeat_subfamily){
  gene <- gene
  repeat_subfamily <- repeat_subfamily
  
  lm_eqn <- function(df){
    m <- lm(log2GeneTPM ~ log2RepeatTPM, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(unname(coef(m)[1]), digits = 2),
                          b = format(unname(coef(m)[2]), digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
  }
  gene_tmp <- gene_expression %>% filter(GeneID == gene) %>% reshape2::melt(id.vars = "GeneID",value.name = "GeneTPM", variable.name = "SampleName" )
  repeat_tmp <- repeat_expression %>% filter(RepeatID == repeat_subfamily) %>% reshape2::melt(id.vars = "RepeatID",value.name = "RepeatTPM", variable.name = "SampleName" )
  combine <- full_join(gene_tmp,repeat_tmp,by="SampleName") %>% 
    mutate(log2RepeatTPM = log2(RepeatTPM),
           log2GeneTPM = log2(GeneTPM )) %>% 
    drop_na()
  r <- with(combine,cor(log2GeneTPM,log2RepeatTPM))
  combine %>% 
    ggplot(aes(x=log2(RepeatTPM),y=log2(GeneTPM)))+
    geom_point(alpha=0.2)+
    geom_smooth(method=lm , color="red", fill="grey", se=TRUE)+
    labs(x=paste0("log2(",repeat_subfamily," TPM)"),y=paste0("log2(",gene," TPM)"),subtitle = paste0("r = ",r))+
    ##  geom_text(x = 1, y = 6, label = lm_eqn(tmp3), size=10,parse = TRUE)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text( size=rel(1.8)),
          axis.title.y = element_text( size=rel(1.8)),
          plot.title = element_text(hjust = 0.5,size=rel(2),face = "bold"),
          plot.subtitle = element_text(hjust = 1,size=rel(1.2)),
          axis.text.x = element_text(angle=0,vjust=0),
          axis.text= element_text(size=rel(1.5)),
          legend.text=element_text(size=rel(1.5)),
          legend.title = element_text(size=rel(1.5)),
          plot.margin = margin(0.3,1,0.3,0.3, "cm"),
          strip.text.x = element_text(size = rel(1.8)),
          strip.placement = "outside",strip.background = element_blank(),
          legend.position = "bottom"
    )+guides(fill=guide_legend(title=""))
}

#### from guangchuang yu github
parse_ratio <- function(ratio) {
  ratio <- sub("^\\s*", "", as.character(ratio))
  ratio <- sub("\\s*$", "", ratio)
  numerator <- as.numeric(sub("/\\d+$", "", ratio))
  denominator <- as.numeric(sub("^\\d+/", "", ratio))
  return(numerator/denominator)
}
get_enrich_count <- function(ratio) {
  ratio <- sub("^\\s*", "", as.character(ratio))
  ratio <- sub("\\s*$", "", ratio)
  denominator <- as.numeric(sub("^\\d+/", "", ratio))
  return(denominator)
}




plot_screen_rank_total <- function(screen_df,genelist,sort_by = "G418.rank",min= -1, max = 1) {
  require(ComplexHeatmap)
  genelist <- str_remove(genelist,"U2AF1")
  plot_screen_heatmap <- function(df,genelist,title="",index = ind){
    # df <- L1_sup
    mat <- df %>% 
      filter(id %in% genelist) %>% 
      select(id,contains("rank")) %>% 
      arrange(id) %>% 
      mutate( across( contains("rank")  , ~ case_when(.x >=10000 ~ " > 10000",
                                                      .x >= 5000 ~ "5000 - 10000",
                                                      .x >= 1000 ~ "1000 - 5000",
                                                      .x >= 500 ~ "500 - 1000",
                                                      .x >= 200 ~ "200 - 500",
                                                      .x >= 100 ~ "100 - 200",
                                                      .x >=50 ~ "50 - 100",
                                                      .x >= 10 ~ "10 - 50",
                                                      .x >= 1 ~ "1 - 10",
                                                      TRUE ~ " > 10000")
      )) %>%
      dft() 
    unique(mat[1,])
    aa <- colnames(mat)
    colors <- structure( c("#de425b","#ea6e79","#f49399","#fbb8ba","#ffdbdc","#b5b9cf","#9197b8","#6d77a1","#49598a") %>% rev()  , names = c(" > 10000",  "5000 - 10000",  "1000 - 5000",  "500 - 1000",  "200 - 500",  "100 - 200",  "50 - 100",  "10 - 50",  "1 - 10")) %>% factor()
    Heatmap(t(mat)[ind,],
            col = colors,
            name = "rank",
            row_names_side = "left",
            column_title = title,
            rect_gp = gpar(col = "white", lwd = 0.5)
    )
  }
  #### sorted by rank
  ind <- screen_df %>% 
    filter(id %in% genelist) %>% 
    arrange(id) %>% 
    pull(!!{sort_by}) %>% 
    order()
  
  ht1 <- plot_screen_heatmap(L1_sup,genelist,"L1 sup",ind)
  ht2 <- plot_screen_heatmap(L1_act,genelist,"L1 act",ind)
  ht3 <- plot_screen_heatmap(ERV_sup,genelist,"LTR sup",ind)
  ht4 <- plot_screen_heatmap(ERV_act,genelist,"LTR act",ind)
  ht_list <- ht1+ht2+ht3+ht4
  
  draw(ht_list, auto_adjust = FALSE)
  
}

read_mutate <- function(DESeq2RepeatMasker,prefix){
  tmp <- read_tsv(file = DESeq2RepeatMasker) %>% select(rowname,log2FoldChange,pvalue,padj)
  tmp$padj[is.na(tmp$padj)] <- 1 
  tmp %>% mutate(sample=rep(prefix,length(rowname)))
}



plot_individual_Topfamily <- function(UniqueDESeqDf,lfc.topn = 10, copynumber.filter = 20,direction = "both"){
  tmp <- UniqueDESeqDf %>% 
    mutate(subfamily = str_split_fixed(rowname,"_chr|_dup",2)[,1]) %>% 
    filter(!str_detect(rowname,"\\(")) %>% 
    group_by(subfamily) %>% 
    mutate( medianlfc = median(log2FoldChange) ) 
  top_n <- tmp %>% 
    summarize( medianlfc = median(log2FoldChange),
               copy = n()) %>% 
    filter(copy > copynumber.filter ) %>% 
    top_bottom(lfc.topn,medianlfc) %>% 
    mutate( name.withnumber = glue("{subfamily} \n ( n = {copy})"))
  plotdf <- left_join(top_n,tmp)
  if (direction == "both") {
      plotdf <- plotdf 
    }else if (direction == "up"){
      plotdf <- plotdf %>% filter(medianlfc >= 0)
    }else if (direction == "down"){
      plotdf <- plotdf %>% filter(medianlfc <= 0)
    }
  increase_number <- top_n %>% arrange((medianlfc)) %>% filter(medianlfc < 0) %>% nrow()
  change <- case_when(direction == "both" ~ "increased/decreased",
                      direction == "up" ~ "increased",
                      direction == "down" ~ "decreased",
                      TRUE ~ direction)
  p <- plotdf %>%
    ungroup() %>% 
    #filter(pvalue < 0.01) %>% 
    mutate(name.withnumber = fct_reorder(name.withnumber,log2FoldChange,median)) %>% 
    ggplot(aes(x=name.withnumber,y = log2FoldChange))+
    geom_violin()+
    geom_boxplot(outlier.colour = "black",width = 0.3)+
    geom_hline(yintercept = 0,color ="red",alpha= 0.3,linetype = 1,size = 1)+
    labs(caption = glue("top {lfc.topn} {change} subfamilies with expressiable copies greater than {copynumber.filter}\n
                        n represents number of experssiable copies in each subfamily"),
         x ="")
  if(direction == "both"){ p <- p+geom_vline(xintercept = increase_number + 0.5,color ="red",alpha= 0.3,linetype = 1,size = 1)}
  return(p)
}

require(tidytext)
plot_individual_subfamily_updown <- function(UniqueDESeqDf , pvalue.type = "pvalue" ,pvalue.cutoff = 0.01,copynumber.filter=10){
  x <- UniqueDESeqDf
  plotdf <- x %>% 
    mutate(subfamily = str_split_fixed(rowname,"_chr|_dup",2)[,1]) %>% 
    filter(!str_detect(rowname,"\\(")) %>% 
    filter( .[[pvalue.type]] < pvalue.cutoff) %>% 
    group_by(subfamily) %>% 
    filter( n() >= copynumber.filter) %>% 
    mutate(updown = ifelse(log2FoldChange > 0 ,"up","down")) %>% 
    group_by(subfamily,updown) %>% 
    summarize(count = n()) %>%
    ungroup %>% 
    complete(subfamily,updown, fill = list(count = 0)) %>% 
    mutate(count = ifelse( updown == "down", -count,count )) 
  
  factor_df <- plotdf %>% 
    ungroup() %>% 
    group_by(subfamily) %>% 
    summarize(total = sum(count)) %>% 
    mutate(subfamily = fct_reorder(subfamily, total,mean))
  
  total_df <- left_join(factor_df,plotdf) %>% 
    mutate(subfamily = factor(subfamily,levels = 	levels(factor_df$subfamily))) %>% 
    ungroup
  #filter(count > copynumber.filter) %>%
  #ungroup() %>% 
  #mutate(subfamily = reorder_within(subfamily,count,updown)) 
  
  
  if(nrow(plotdf) != 0){
    total_df %>% 
      ggplot(aes(x=count,y=subfamily,fill = updown)) +
      geom_bar(stat="identity")+
      geom_point(aes(x=total,y = subfamily),shape = 18,size = 3)+
      scale_x_continuous(labels = abs)+
      ggsci::scale_fill_aaas()+
      xlab("The number of instances")+ ylab("")+
      labs(caption = glue("copynumber filter : {copynumber.filter} \n {pvalue.type} < {pvalue.cutoff}"))
    
  }else { ggplot() }
}


RepeatMAplot_family <- function(DESeqdf,p="pvalue",pcutoff=0.05,family_name="L1",filterExpression=0,title="",label_top_n=5){
  df <- as.data.frame(DESeqdf)
  df$padj[is.na(df$padj)] <- 1
  colnames(df)[1] <- "symbol"
  family <- read_tsv("/LiuLab/reference/Human/GRCh38/TE/Repeat_Masker_hg38.class",col_names=F)
  colnames(family) <- c("subfamily","class","family")
  repeat_family_name <- family$family %>% unique
  family_list <- family %>% filter(family == family_name) %>% pull(subfamily)
    aa <- ifelse(df$symbol %in% family_list & df[,colnames(df)==p]< pcutoff,paste0("significant ",family_name),ifelse(df$symbol %in% family_list  ,family_name,paste0("others")))
    aa[is.na(aa)] <- paste0("others")
    data <- df %>% add_column(group = aa)
    data$group <- factor(data$group,levels=c(paste0("others"),family_name,paste0("significant ",family_name)))
    ### label data frame
    labeldf <- get_top_label(data=data,group=paste0("significant ",family_name),label_top_n)
    
    plot <- ggplot(data=data,mapping = aes(x=baseMean,y=log2FoldChange,color=group,size=group,alpha=group))+
      geom_point()+
      geom_hline(yintercept = 0.0, color = "red", size = 1.5,alpha=0.3)+
      scale_colour_manual(values = c("grey","blue", "red"))+
      scale_alpha_manual(values = c(0.1,0.7,1))+
      scale_size_manual(values = c(0.5,1.5,2.5))+
      ggrepel::geom_text_repel(data=labeldf,label=labeldf$symbol,size=7, force = 2,segment.size=0.2,segment.alpha=0.95,segment.color="black",color="black")+
      scale_x_continuous(labels = comma,trans = "log10")+
      #ylim(-1.5,1.5)+
      xlab("RNA-seq reads (normalized counts)")+
      ylab("Gene expression change (log2[Treatment/control])") +
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
  )+
      guides(color=guide_legend(title="",override.aes = list(size=5)),alpha="none",size="none")+
      ggtitle(paste0(title," , ",p," < ",pcutoff,", ",family_name))
  
  #p2 <- ggMarginal(p,type = "density",margins = "y",groupColour  = T,groupFill = T)
 
}

foldchange_plot <- function(df){
  ggplot(df,aes(x=LINE1, y=sample, size=padj, fill=log2FoldChange)) +
    geom_point(alpha=1, shape=21, color="black") +
    scale_size(range = c(8, .1), name="padj",breaks =c(0.001,0.01,0.1,0.5,1)) +
    scale_fill_gradient2( low = "darkblue",
                          mid = "white",
                          high = muted("red"),
                          breaks=c(-4,-2,-1,-0.5,0,0.5,1,1.5,2,3,4),
                          )+
    # theme_minimal() +
    theme(#panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        axis.title.x = element_text( size=rel(1)),
        axis.title.y = element_text( size=rel(1)),
        plot.title = element_text(hjust = 0.5,size=rel(1),face = "bold"),
        plot.subtitle = element_text(hjust = 1,size=rel(1)),
	      axis.text.x = element_text(angle=90,vjust=0.5),
        axis.text= element_text(size=rel(1)),
        legend.text=element_text(size=rel(1)),
        legend.title = element_text(size=rel(1)),
 	      plot.margin = margin(0.3,1,0.3,0.3, "cm"),
 	      strip.text.x = element_text(size = rel(1)),
		    strip.placement = "outside",strip.background = element_blank(),
        legend.position = "right"
  )+
    xlab("")+ylab("")
}
