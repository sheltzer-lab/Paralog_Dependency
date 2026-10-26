################################
########### Figure 4 ###########
################################
# Project: CYCLOPS are singleton 
# Author: Yi

##################################
########### Figure 4F ############
##################################
# library
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(ggvenn)

# load the data
sig_df = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/Paralog dependency/codes/Paralogs identity/dep_summary_aneuploidy_chr_loss.xlsx")

# sig CYCLOPS target from my analysis
# the chormosome loss events 
sig_df %>% 
  dplyr::filter(status == "CYCLOPS") %>% 
  dplyr::select(dep_gene,aneuploid_loss_chr) -> CYCLOPS_df


# GSEA analysis from our GSEA analysis 
# function analysis 
library(clusterProfiler)

GO_analysis = clusterProfiler::enrichGO(CYCLOPS_df$dep_gene,
                                        'org.Hs.eg.db', ont="BP" ,
                                        keyType = "SYMBOL",pvalueCutoff = 0.05)

View(GO_analysis[1:20,])
GO_analysis_df = as.data.frame(GO_analysis)
writexl::write_xlsx(GO_analysis_df, "GO_CYCLOPS.xlsx")

# plot the GSEA analysis
library(MetBrewer)
GO_analysis_df %>% 
  .[c(1,2,3,5,6,7,8,10,12,13),] %>% 
  ggplot(aes(x = reorder(Description, -qvalue),y = -log10( qvalue))) +
  geom_point(aes(size = Count ,color = qvalue)) +
  coord_flip() +
  geom_rangeframe() +
  theme_tufte() +
  scale_color_gradient2(high ="#083D77", mid ="#ebebd3" ,low = "#E2725B",midpoint = 0.018)

# library
library(tidyverse)
library(ggthemes)
library(ggpubr)

# colors
colors = c("#EF8536","#3A76AF")
# dataset
# load hg38 genomes
library(biomaRt)

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

mart <- useMart(biomart="ensembl",
                dataset="hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position",
                "end_position", "strand", "gene_biotype", "description","band")

# Fetch the gene information
gene_info <- getBM(attributes = attributes, mart = mart,useCache = TRUE)

library(MetBrewer)
# check the gene annotation and distribution
gene_info %>% 
  dplyr::filter(chromosome_name %in% c(1:22,"Y","X")) %>% 
  dplyr::filter(hgnc_symbol != "") %>% 
  distinct(hgnc_symbol, .keep_all = T) %>% 
  group_by(gene_biotype) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>% 
  dplyr::mutate(gene_biotype = factor(gene_biotype, levels = gene_biotype))->p
p %>% 
  ggplot(aes(x =gene_biotype, y = count)) +
  geom_histogram(stat = "identity", width = 0.7,aes(color = gene_biotype, fill = gene_biotype), alpha = 0.8) +
  geom_rangeframe() +
  theme_tufte()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values  = met.brewer("Hokusai3", n = 33 ,direction = 1)) +
  scale_color_manual(values = met.brewer("Hokusai3" ,n = 33,direction = 1))

# include all the genes from the 
# check how many genes with paralogs
# focus on the with protein-coding
gene_info %>% 
  dplyr::filter(chromosome_name %in% c(1:22,"Y","X")) %>% 
  dplyr::filter(hgnc_symbol != "") %>% 
  distinct(hgnc_symbol, .keep_all = T) %>% 
  dplyr::filter(gene_biotype == "protein_coding")-> all_chr_gene


# gene with parlogus 
paralogs <- getBM(attributes = c("hsapiens_paralog_associated_gene_name",
                                 "ensembl_gene_id",
                                 "hsapiens_paralog_ensembl_gene",
                                 "hsapiens_paralog_perc_id"),
                  filters = "ensembl_gene_id", 
                  values =all_chr_gene$ensembl_gene_id, 
                  mart = mart)


# integrate the genes_names with paralogs data
# filter out the no identical score pair
# get the gene with pairs
all_chr_gene %>% 
  left_join(paralogs, by = "ensembl_gene_id") %>% 
  dplyr::filter(!is.na(hsapiens_paralog_perc_id)) %>% 
  dplyr::filter(hsapiens_paralog_associated_gene_name != "") %>% 
  dplyr::mutate(para_pair_for = paste(hgnc_symbol,hsapiens_paralog_associated_gene_name,sep = "_" )) %>% 
  dplyr::mutate(para_pair_rev = paste(hsapiens_paralog_associated_gene_name,hgnc_symbol,sep = "_" ))-> paralog_gene_pair

# check the head of paralogs pair
View(head(paralog_gene_pair))

# there is an issue that the paralogs have the forward and reverse pair
# this goona cause the duplicate he count of genes
# we need to avoid the duplicate count of paralogs 
# identical score is different calculate in the same pair of paralogs 
# we goona gonna calculate the mean of identical score

for_para = paralog_gene_pair$para_pair_for
rev_para = paralog_gene_pair$para_pair_rev

# read the paralog pair both side
reformed_para_cmb = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/Paralog dependency/codes/Paralogs identity/paralog_pair_both_side.csv")
reformed_para_cmb = reformed_para_cmb[,-1]

## To get the unique paralog pairs
# remove the duplicates
# index_rev corresponding on paralog_gene_pair
reformed_para_cmb %>% 
  dplyr::filter(index_for < index_rev | is.na(index_rev) ) %>% 
  distinct(para_name, .keep_all = T)-> reformed_para_cmb_flt # distinct the duplication # no duplicates paralog pair

# define the paralgs by identical score
# calcluate the famaily size
reformed_para_cmb %>% # paralog gene family
  distinct(para_name, .keep_all = T) %>% 
  dplyr::filter(mean_identical_score >=20) %>% 
  group_by(hgnc_symbol) %>% 
  summarise(family_size = n()) -> family_size_sum

# filter the paralogs
reformed_para_cmb_flt %>% 
  dplyr::filter(mean_identical_score >= 20) %>% 
  dplyr::filter(gene_biotype == "protein_coding") %>% 
  left_join(family_size_sum, by = "hgnc_symbol") %>% 
  dplyr::filter(family_size < 20) -> paralog_pair_flt


# get the CYClOPS gene
# read the significant target 
sig_target = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/Paralog dependency/codes/Paralogs identity/dep_summary_aneuploidy_chr_loss.xlsx")

sig_target %>% 
  dplyr::filter(status == "CYCLOPS") -> df_CYCLOPS

#
paralog_num = c()
for (gene in df_CYCLOPS$dep_gene) {
  if(!(gene %in% family_size_sum$hgnc_symbol)){
    family_size = 0
    # add the family size
    paralog_num = c(paralog_num, family_size)
  }else{
    family_size_sum %>% 
      dplyr::filter(hgnc_symbol == gene) %>% 
      distinct(hgnc_symbol, .keep_all = T)-> tpm_df
    family_size = tpm_df$family_size
    # add the family size
    paralog_num = c(paralog_num, family_size)
    
    
  }
}

CYCLOPS_family_size = data.frame(df_CYCLOPS$dep_gene, paralog_num)

CYCLOPS_family_size %>% 
  dplyr::mutate(paralog_num = paralog_num +1) -> CYCLOPS_family_size_mdf



################ non-sig CYCLOPS #############
####### gene on the chromosome loss arm ##########
all_chr_gene %>% 
  dplyr::filter(!(chromosome_name %in% c("Y", "X"))) %>% 
  dplyr::mutate(chr_band = paste(chromosome_name, str_sub(band,1,1), sep = "")) %>% 
  dplyr::filter(chr_band %in% c("1p", "3p", "4p", "4q", "5q", "6p", "6q", "8p", "9p", "9q",
                                "10p", "10q", "11p", "11q", "12p", "13q", "14q", "15q", "16p",
                                "16q", "17p", "18p", "18q", "19p", "19q", "20p", "21q", "22q")) %>% 
  dplyr::select(hgnc_symbol) %>% 
  unlist(use.names = F) -> loss_chr_gene_lists

# non_sig_cyclops
non_sig_cyclops = loss_chr_gene_lists[-which(loss_chr_gene_lists %in% CYCLOPS_family_size$df_CYCLOPS.dep_gene)]

#family_size_sum
non_sig_cyclops_df = data.frame(non_sig_cyclops)
colnames(non_sig_cyclops_df) = "hgnc_symbol"
non_sig_cyclops_df %>% 
  left_join(family_size_sum,  by = "hgnc_symbol") %>% 
  dplyr::mutate(family_size_new = ifelse(is.na(family_size), 0, family_size)) %>% 
  dplyr::filter(family_size_new <=20) -> non_sig_cyclops_gene_family_size

non_sig_cyclops_gene_family_size %>% 
  dplyr::mutate(family_size_new = family_size_new + 1) %>% 
  dplyr::select(hgnc_symbol, family_size_new) -> non_sig_cyclops_gene_family_size_mdf


# compare together to check the family size between sig and non-sig CYCLOPS based on the famliy size 
colnames(CYCLOPS_family_size_mdf) = c("hgnc_symbol", "Family_size")
CYCLOPS_family_size_mdf$condition = "Significant_CYCLOPS"
colnames(non_sig_cyclops_gene_family_size_mdf) = c("hgnc_symbol", "Family_size")
non_sig_cyclops_gene_family_size_mdf$condition = "Non-Sig_CYCLOPS"

# combine 2 tables 
CYCLOPS_cmb = rbind(CYCLOPS_family_size_mdf, non_sig_cyclops_gene_family_size_mdf)



plot_df <- CYCLOPS_cmb %>%
  dplyr::count(condition, Family_size, name = "size_count") %>%
  dplyr::mutate(
    condition = factor(
      condition,
      levels = c("Significant_CYCLOPS", "Non-Sig_CYCLOPS")
    )
  ) %>%
  tidyr::complete(
    condition,
    Family_size,
    fill = list(size_count = 0)
  ) %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(
    total = sum(size_count),
    freq = size_count / total
  ) %>%
  dplyr::ungroup()

#############################################
################# Figure 4G #################
#############################################
# Plot the figure
ggplot(plot_df, aes(x = factor(Family_size), y = freq)) +
  
  geom_histogram(aes( fill = condition, color = condition), 
                 stat = "identity",
                 position = position_dodge(width = 0.9),
                 alpha= 0.8,width = 0.6) +
  geom_rangeframe() + 
  theme_tufte() +
  geom_rangeframe() +
  theme_tufte() +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_discrete(breaks = as.character(seq(1, 21, by = 2)), drop = FALSE) +
  labs(x = "Family size", y = "Frequency")

# chisqr test
plot_df %>% 
  dplyr::filter(Family_size == 1) -> df_plot_family_size1

#p = 0.001263
chisq.test(
  matrix(
    c(18, 10,
      3770, 7433),
    nrow = 2,
    byrow = TRUE
  )
) 

#############################################
################# Figure 4H #################
#############################################
###################  Dependency difference between 2 groups ####################
# dependency filter 
# To judge that paralog could create dependency, we integrate the Depmap data and expression data from CCLE

# cell line info
metasheet_cellline = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /Model.csv")

# filter the cell line are aneuploid
aneuploid_dep = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/aneuploidy_scores.csv")

aneuploid_dep %>% 
  dplyr::filter(Aneuploidy.score >= 7) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> aneuploid_cellline

df_aneu_score = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/arm_call_scores.csv")


# Gene dependency 24Q2 data set
df_gene_dep = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /CRISPRGeneEffect.csv")

# change the name of data
colnames(df_gene_dep)[1] = "DepMap_ID"
dim(df_gene_dep)

# gene the gene name 
dim(df_gene_dep)
genelist_Dep = sub("\\..*","",colnames(df_gene_dep)[2:dim(df_gene_dep)[2]])

# replace gene to colnames
colnames(df_gene_dep)[2:dim(df_gene_dep)[2]] = genelist_Dep
head(df_gene_dep[,1:6])

df_gene_dep_flt = df_gene_dep[df_gene_dep$DepMap_ID %in% aneuploid_cellline, ]
Gene_dependency_score = function(x){
  # get the gene
  gene = x
  # get the dependency
  gene_dep = df_gene_dep_flt[,which(colnames(df_gene_dep_flt) == gene)]
  
  if(!is_empty(gene_dep)){
    # calculate the median of the dependency
    median_dep = median(gene_dep[!is.na(gene_dep)])
    
    
  }else{
    median_dep = NA
  }
  
  return(median_dep)
  
}

# get the dependency
#dependency
sig_CYCLOPS_dep = lapply(CYCLOPS_family_size$df_CYCLOPS.dep_gene, Gene_dependency_score)
sig_CYCLOPS_dep_cmb = do.call(rbind, sig_CYCLOPS_dep)
sig_CYCLOPS_dep_cmb = as.data.frame(sig_CYCLOPS_dep_cmb)
sig_CYCLOPS_dep_cmb$condition = "sig"
colnames(sig_CYCLOPS_dep_cmb)[1] = "count"

# non_sig
nonsig_singleton_dep = lapply(non_sig_cyclops, Gene_dependency_score)
nonsig_singleton_dep_cmb = do.call(rbind, nonsig_singleton_dep)
nonsig_singleton_dep_cmb = as.data.frame(nonsig_singleton_dep_cmb)
nonsig_singleton_dep_cmb$condition = "nonsig"
colnames(nonsig_singleton_dep_cmb)[1] = "count"


# combine the file
dep_cmb_df = rbind(sig_CYCLOPS_dep_cmb, nonsig_singleton_dep_cmb)


dep_cmb_df %>% 
  dplyr::mutate(condition = factor(condition,levels = c("sig", "nonsig"))) %>% 
  ggplot(aes(x = condition,y = count)) +
  geom_boxplot(
    aes(color = condition, fill = condition),
    width = 0.25,
    outlier.shape = NA, outliers = F,
    alpha = 0.4
  ) +
  geom_rangeframe() +
  theme_tufte() +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(x = "", y = "dependency score") +
  ylim(c(-2,1))

#
#p-value = 9.78e-15
wilcox.test(dep_cmb_df[dep_cmb_df$condition == "sig",]$count, 
            dep_cmb_df[dep_cmb_df$condition == "nonsig",]$count)

#############################################
################# Figure 4I #################
#############################################
# CYCLOPS expression and dependency correlation
sig_cyclops_dep_exp_cor = function(x){
  
  # input list 
  df_tpm = cyclops_loss_chr_gene_df[x,]
  
  # input  gene_dep
  gene_dep = df_tpm$dep_gene
  
  # input _loss gene
  gene_loss = df_tpm$dep_gene
  
  # get the expresssion of gene_loss
  exp_flt_df = CCLE_exp[,c(1,which(colnames(CCLE_exp) == gene_loss))]
  
  # get the dependency data
  dep_flt_list = df_gene_dep[,c(1, which((colnames(df_gene_dep) ==gene_dep )))]
  
  
  exp_flt_df %>% 
    left_join(dep_flt_list, by = "DepMap_ID") %>% 
    dplyr::filter(.[[2]] > 1) %>% 
    dplyr::filter(DepMap_ID %in% aneuploid_cellline) -> df_flt_exp_dep
  
  
  
  cor_score = cor.test(df_flt_exp_dep[,2],df_flt_exp_dep[,3])$estimate
  df_return =data.frame(gene_dep,gene_loss ,cor_score)
  return(df_return)
} 

# apply for all the sig paralogs 
cor_exp_dep_list = lapply(c(1:dim(cyclops_loss_chr_gene_df)[1]), sig_cyclops_dep_exp_cor)

cor_exp_dep_list_cmb = do.call(rbind,cor_exp_dep_list)

cor_exp_dep_list_cmb %>% 
  ggplot(aes(x = cor_score, y = ..scaled..)) +
  geom_density(fill = "#5271A8", color = "#5271A8", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray49") +
  #geom_rangeframe() +
  theme_tufte() +
  xlim(c(-0.5,1)) +
  labs(x = "Correlation", y = "Scaled density") 



#############################################
################# Figure 4J #################
#############################################
############## Check for the expression ##################

# retrieve gene expression data
CCLE_exp = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")

# change the name of CCL2_exp data
colnames(CCLE_exp)[1] = "DepMap_ID"
dim(df_gene_dep)

# gene the gene name 
dim(df_gene_dep)
genelist_exp = sub("\\..*","",colnames(CCLE_exp)[2:dim(CCLE_exp)[2]])

# replace gene to colnames
colnames(CCLE_exp)[2:dim(CCLE_exp)[2]] = genelist_exp

###########################################################
### sig CYCLOPS gene expression in aneuploid cancer cell###

SIG_CYCLOPS_exp = function(x){
  # get the gene name and chromsome postion
  
  df_tmp = sig_CYCLOPS_mdf[x,]
  
  # gene  
  
  gene_name = df_tmp[,1]
  
  # chromosome position
  
  chr_pos = df_tmp[,4]
  
  # get the certainc chromosome info
  loss_chr_condition = aneuploid_status_flt[,c(1, which(colnames(aneuploid_status_flt) == chr_pos))]
  #look a the certain chrosmome loss cell
  chr_loss_cells = loss_chr_condition[loss_chr_condition[,2] == -1, ]$X
  
  #look a the certain chrosmome neutral cell
  chr_neutral_cells = loss_chr_condition[loss_chr_condition[,2] == 0, ]$X
  
  # get the expression level
  gene_exp_flt = CCLE_exp[,c(1, which(colnames(CCLE_exp) ==gene_name ))]
  
  # get the loss expression 
  gene_exp_flt %>% 
    dplyr::filter(DepMap_ID %in% chr_loss_cells) %>% 
    .[,2] %>% unlist(use.names = F) -> loss_exp
  
  # get the neutral expression 
  gene_exp_flt %>% 
    dplyr::filter(DepMap_ID %in% chr_neutral_cells) %>% 
    .[,2] %>% unlist(use.names = F) -> neutral_exp
  
  # normalised the expression
  median_exp_loss = median(loss_exp[!is.na(loss_exp)])
  
  median_exp_neutral = median(neutral_exp[!is.na(neutral_exp)])
  
  normm_loss_exp= median_exp_loss- median_exp_neutral
  
  norm_neutral_exp = median_exp_neutral - median_exp_neutral
  
  df_return = data.frame(gene_name,normm_loss_exp, norm_neutral_exp)
  
  return(df_return)
  
}

sig_CYCLOPS_exp_norm = lapply(c(1:dim(sig_CYCLOPS_mdf)[1]), SIG_CYCLOPS_exp)  
sig_CYCLOPS_exp_norm_cmb = do.call(rbind, sig_CYCLOPS_exp_norm)


sig_CYCLOPS_exp_norm_cmb %>% 
  gather(normm_loss_exp:norm_neutral_exp, key = "condition", value = "value") %>% 
  dplyr::mutate(condition = factor(condition, levels = c("normm_loss_exp", "norm_neutral_exp"))) -> sig_CYCLOPS_exp_mdf

sig_CYCLOPS_exp_mdf %>% 
  ggplot(aes(x = condition,y = value)) +
  geom_boxplot(
    aes(color = condition, fill = condition),
    width = 0.25,
    outlier.shape = NA, outliers = F,
    alpha = 0.4
  ) +
  geom_rangeframe() +
  theme_tufte() +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  scale_y_continuous(breaks = seq(-1, 0.75, by = 0.25)) +
  coord_cartesian(ylim = c(-1, 0.75)) +
  labs(x = "", y = "Relative expression") -> P3



################### non sig #########################
NONSIG_CYCLOPS_exp = function(x){
  # get the gene name and chromsome postion
  
  df_tmp = nonsig_CYCLOPS_mdf[x,]
  
  # gene  
  
  gene_name = df_tmp[,1]
  # print(gene_name)
  
  # chromosome position
  
  chr_pos = df_tmp[,4]
  
  # get the certainc chromosome info
  loss_chr_condition = aneuploid_status_flt[,c(1, which(colnames(aneuploid_status_flt) == chr_pos))]
  #look a the certain chrosmome loss cell
  chr_loss_cells = loss_chr_condition[loss_chr_condition[,2] == -1, ]$X
  
  #look a the certain chrosmome neutral cell
  chr_neutral_cells = loss_chr_condition[loss_chr_condition[,2] == 0, ]$X
  
  # get the exp score
  gene_exp_flt = CCLE_exp[,c(1, which(colnames(CCLE_exp) ==gene_name ))]
  if(!is.null(dim(gene_exp_flt))){
    # get the loss exp 
    gene_exp_flt %>% 
      dplyr::filter(DepMap_ID %in% chr_loss_cells) %>% 
      .[,2] %>% unlist(use.names = F) -> loss_exp
    
    # get the neutral exp 
    gene_exp_flt %>% 
      dplyr::filter(DepMap_ID %in% chr_neutral_cells) %>% 
      .[,2] %>% unlist(use.names = F) -> neutral_exp
    
    # normalised the exp
    median_exp_loss = median(loss_exp[!is.na(loss_exp)])
    
    median_exp_neutral = median(neutral_exp[!is.na(neutral_exp)])
    
    normm_loss_exp = median_exp_loss- median_exp_neutral
    
    norm_neutral_exp = median_exp_neutral- median_exp_neutral
  }else{
    normm_loss_exp = NA
    
    norm_neutral_exp = NA
  }
  
  
  df_return = data.frame(gene_name,normm_loss_exp, norm_neutral_exp)
  
  return(df_return)
  
}

nonsig_CYCLOPS_exp_norm = lapply(c(1:dim(nonsig_CYCLOPS_mdf)[1]), NONSIG_CYCLOPS_exp)  
nonsig_CYCLOPS_exp_norm_cmb = do.call(rbind, nonsig_CYCLOPS_exp_norm)

nonsig_CYCLOPS_exp_norm_cmb %>% 
  gather(normm_loss_exp:norm_neutral_exp, key = "condition", value = "value") %>% 
  dplyr::mutate(condition = factor(condition, levels = c("normm_loss_exp", "norm_neutral_exp"))) -> non_sig_CYCLOPS_exp_mdf
non_sig_CYCLOPS_exp_mdf %>% 
  ggplot(aes(x = condition,y = value)) +
  geom_boxplot(
    aes(color = condition, fill = condition),
    width = 0.25,
    outlier.shape = NA, outliers = F,
    alpha = 0.4
  ) +
  geom_rangeframe() +
  theme_tufte() +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  scale_y_continuous(breaks = seq(-1, 0.75, by = 0.25)) +
  coord_cartesian(ylim = c(-1, 0.75)) +
  labs(x = "", y = "Relative expresssion") -> P4

# Statistic
# p= 4.905e-08
wilcox.test(non_sig_CYCLOPS_exp_mdf[non_sig_CYCLOPS_exp_mdf$condition == "normm_loss_exp", ]$value, 
            sig_CYCLOPS_exp_mdf[sig_CYCLOPS_exp_mdf$condition == "normm_loss_exp", ]$value)


ggarrange(P3,P4)


# Project: 72 significant interactors and correlation
# Date: Nov.03.2025
# Author: Yi 


# library
library(tidyverse)
library(ggthemes)
library(ggpubr)
#library(biomaRt)

# load data
# dataset
# CCLE aneuploid status annotation
aneuploid_dep = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/aneuploidy_scores.csv")
# CCLE cell line with each chromsome arm gain or loss info
df_aneu_score = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/arm_call_scores.csv")
#cell line annotation
cell_line_id = read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/Cell_lines_annotations_20181226.txt")

# get aneuploid cells
# filter the cell line are aneuploid
aneuploid_dep %>% 
  dplyr::filter(Aneuploidy.score >= 7) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> aneuploid_cellline

# get the significant partner
# read the file all significant 37 significant paralogs
sig_dep_df = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/Paralog dependency/codes/Paralogs identity/dep_summary_aneuploidy_chr_loss.xlsx")


# get the cell interactors dependency
sig_dep_df %>% 
  dplyr::filter(status == "interaction_dep") -> int_df_flt

# load hg38 genomes
library(biomaRt)
mart <- useMart(biomart="ensembl",
                dataset="hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position",
                "end_position", "strand", "gene_biotype", "description","band")

# Fetch the gene information
gene_info <- getBM(attributes = attributes, mart = mart,useCache = TRUE)
# filter the chromosome 
gene_info %>% 
  dplyr::filter(chromosome_name %in% c(1:22,"Y","X")) %>% 
  dplyr::filter(hgnc_symbol != "") %>% 
  distinct(hgnc_symbol, .keep_all = T) %>% 
  dplyr::filter(gene_biotype == "protein_coding") %>% 
  dplyr::mutate(chr_pos = ifelse(chromosome_name %in% c("X", "Y"), paste("chr", chromosome_name, sep = ""),
                                 paste(paste("chr",chromosome_name, sep = "" ), str_sub(band, 1,1),sep = "_"))) %>% 
  dplyr::select(hgnc_symbol, chr_pos)-> all_chr_gene
# library
library(rbioapi)

# write function to get the gene interactors
interactor_number_likelihood = function(x){
  # input the data 
  input_df =  int_df_flt[x,]
  
  # dep_gene
  dep_gene = input_df$dep_gene
  
  # loss chromosome
  loss_chr_spec = input_df$aneuploid_loss_chr
  
  # get the dep gene 
  proteins_mapped_1 <- rba_string_map_ids(ids = dep_gene,
                                          species = 9606)
  
  
  proteins_mapped_query = proteins_mapped_1$stringId
  
  int_partners <- rba_string_interaction_partners(ids = proteins_mapped_query,
                                                  species = 9606,
                                                  required_score = 500,
                                                  limit = 200) # set the limit interactor number 200
  # cutoff at 0.7 and assign the chromosome position
  int_partners %>% 
    dplyr::filter(score > 0.5) %>% 
    dplyr::rename( hgnc_symbol = preferredName_B) %>% 
    left_join(all_chr_gene, by = "hgnc_symbol") -> df_int_partner
  
  # get the summary of number of interactors
  interactor_count = dim(df_int_partner)[1]
  df_int_partner %>% 
    dplyr::filter(chr_pos == loss_chr_spec) -> inter_loss_chr_flt
  
  # get the count the interactors 
  if(!is.null(dim(inter_loss_chr_flt))){
    interactor_loss_count = dim(inter_loss_chr_flt)[1]
  }else{
    interactor_loss_count = 0
  }
  
  
  return_df = data.frame(dep_gene, interactor_count, interactor_loss_count)
  colnames(return_df) = c("gene", "total_inter_count", "loss_chr_inter_count")
  
  return(return_df)
  
  
}

# lapply the function to the data
summary_inter_73 = lapply(c(1:dim(int_df_flt)[1]), interactor_number_likelihood)

summary_inter_73_cmb = do.call(rbind, summary_inter_73)

summary_inter_73_cmb %>% 
  dplyr::mutate(inter_count = loss_chr_inter_count) %>% 
  group_by(inter_count) %>% 
  dplyr::mutate(total_inter = sum(total_inter_count)) %>% 
  dplyr::mutate(specific_inter = sum(loss_chr_inter_count)) %>% 
  ungroup() %>% 
  distinct(inter_count, .keep_all = T) %>% 
  dplyr::mutate(freq_inter = specific_inter/total_inter )  -> df_tpm




###################################################
#################### Figure 4P ####################
###################################################
# WE LOOK AT HYPK 

# get the dep gene 
proteins_mapped_HYPK <- rba_string_map_ids(ids = "HYPK",
                                           species = 9606)


HYPK_mapped_query = proteins_mapped_HYPK$stringId

HYPK_int_partners <- rba_string_interaction_partners(ids = HYPK_mapped_query,
                                                     species = 9606,
                                                     required_score = 500,
                                                     limit = 100) 
HYPK_int_partners$preferredName_B

# library for ggraph
library(igraph)
library(ggraph)
library(graphlayouts)
library(ggforce)
library(scatterpie)

# read HYPK interaction network
HYPK_network = read.delim("HYPK_string_interactions.tsv")
colnames(all_chr_gene)[1] = "X.node1"

HYPK_network %>% 
  left_join(all_chr_gene, by = "X.node1") %>% 
  dplyr::select(X.node1, node2, chr_pos) -> HYPK_network_edt

unique(HYPK_network_edt$chr_pos)
# 
g <- graph_from_data_frame(
  d = HYPK_network_edt[,1:2],      # edges
  directed = FALSE       # or TRUE if directional
)

g <- igraph::simplify(g)

HYPK_network_edt %>% 
  distinct(X.node1, .keep_all = T) -> network_dist

V(g)$grp <- network_dist$chr_pos
V(g)$grp <- factor(V(g)$grp, levels = unique(network_dist$chr_pos))
library(oaqc)
bb <- layout_as_backbone(g, keep = 0.4)
E(g)$col <- F
E(g)$col[bb$backbone] <- T

#
library(MetBrewer)

nodes_to_label <- c("HYPK", "NAA15", "NAA11")

# plot
ggraph(g,
       layout = "manual",
       x = bb$xy[, 1],
       y = bb$xy[, 2]) +
  geom_edge_link0(aes(col = col), width = 0.2) +
  geom_node_point(aes(fill = grp), shape = 21, size = 3) +
  geom_node_text(
    aes(label = ifelse(name %in% nodes_to_label, name, "")),
    repel = TRUE,
    size = 3
  ) +
  geom_mark_hull(
    aes(x, y, group = grp, fill = grp, label=grp),
    concavity = 4,
    expand = unit(2, "mm"),
    alpha = 0.25
  ) +
  
  
  scale_color_manual(values  = c(rep("#0000004D", time = 2), "#e58068","#0000004D" , "#6b7cb9", rep("#0000004D", time = 10))) +
  scale_fill_manual(values  = c(rep("#0000004D", time = 2), "#e58068", "#0000004D", "#6b7cb9", rep("#0000004D", time = 10))) +
  scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
  theme_graph()+
  theme(legend.position = "none")



##################################
########### Figure 4Q ############
##################################

#weigth
library(MetBrewer)
df_tpm$pc = predict(prcomp(~inter_count+freq_inter, df_tpm))[,1]
df_tpm %>% 
  ggplot(aes(x = inter_count,y = freq_inter)) +
  geom_point(aes(color = pc,,size = freq_inter),alpha = 0.8) +
  geom_smooth(method = "lm",alpha = 0.1) +
  geom_rangeframe() +
  theme_tufte() + 
  scale_color_gradient2(low = "#FCAE1E" ,mid = "#E2DDFC",  high = "#543D7B") + 
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(1, 16, by = 4)) +
  labs(x = "number of interator partner on one chromosome", y = "possibility of significant dependent targets")

#8.476e-06
#  cor 
#0.9345387 
cor_result = cor.test(df_tpm$inter_count, df_tpm$freq_inter)



##################################
########### Figure 4R ############
##################################

# project: GSEA analysis of significant interaction targets 

# library
library(tidyverse)
library(ggthemes)
library(clusterProfiler)
library(ggpubr)

# data
# load the data
sig_df = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/Paralog dependency/codes/Paralogs identity/dep_summary_aneuploidy_chr_loss.xlsx")

# retrieve the significant interation genes
sig_df %>% 
  dplyr::filter(status == "interaction_dep") %>% 
  dplyr::select(dep_gene) %>% 
  unlist(use.names = F) ->sig_int_genes

# GSEA analysis from our GSEA analysis 
# function analysis 

GO_analysis = clusterProfiler::enrichGO(sig_int_genes,
                                        'org.Hs.eg.db', ont="BP" ,
                                        keyType = "SYMBOL",pvalueCutoff = 0.05)

View(GO_analysis[1:20,])
GO_analysis_df = as.data.frame(GO_analysis)
writexl::write_xlsx(GO_analysis_df, "Interactor_GO.xlsx")

# add the order ID and easy to pick up the annotaion
GO_analysis_df %>% 
  dplyr::mutate(ID_order = c(1:dim(GO_analysis_df)[1])) -> GO_item_order


GO_item_order %>% 
  .[c(1,2,3,4,5,6,7,9,27,38),] %>% 
  ggplot(aes(x = reorder(Description, -qvalue),y = -log10( qvalue))) +
  geom_point(aes(size = Count ,color = qvalue)) +
  coord_flip() +
  geom_rangeframe() +
  theme_tufte() +
  scale_color_gradient2(high ="#083D77", mid ="#ebebd3" ,low = "#E2725B",midpoint = 0.01)




