#################################
########### Figure 2  ###########
#################################
# Author: Yi
# Date: Oct_2026

# Library
library(tidyverse)
library(ggthemes)
library(ggpubr)

# load hg38 genomes
library(biomaRt)
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
  filter(!is.na(hsapiens_paralog_perc_id)) %>% 
  filter(hsapiens_paralog_associated_gene_name != "") %>% 
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
reformed_para_cmb = read.csv("paralog_pair_both_side.csv")
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


# In this analysis we only focus on the cell line are aneuploidy 
# load aneuploid score dataset

# CCLE aneuploid status annotation
aneuaploid_dep = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/aneuploidy_scores.csv")

# filter the cell line are aneuploid
aneuaploid_dep %>% 
  dplyr::filter(Aneuploidy.score >= 7) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> aneuploid_cellline

# get the aneuploidy condition
df_aneu_score = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/arm_call_scores.csv")

# filter out the naueploidy cells
df_aneu_score %>% 
  dplyr::filter(X %in% aneuploid_cellline) -> df_aneu_score_flt

# dependency filter 
# To judge that paralog could create dependency, we integrate the Depmap data and expression data from CCLE
# Gene dependency 24Q2 data set
df_gene_dep = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /CRISPRGeneEffect.csv")

# cell line info
metasheet_cellline = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /Model.csv")

# change the name of data
colnames(df_gene_dep)[1] = "DepMap_ID"
dim(df_gene_dep)

# gene the gene name 
dim(df_gene_dep)
genelist_Dep = sub("\\.\\..*","",colnames(df_gene_dep)[2:dim(df_gene_dep)[2]])

# replace gene to colnames
colnames(df_gene_dep)[2:dim(df_gene_dep)[2]] = genelist_Dep
head(df_gene_dep[,1:6])

## analysis here: correlation and dep and exp
# retrieve gene expression data from CCLE
CCLE_exp = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")

# change the name of CCL2_exp data
colnames(CCLE_exp)[1] = "DepMap_ID"
#dim(df_gene_dep)

# gene the gene name 
#dim(df_gene_dep)
genelist_exp = sub("\\.\\..*","",colnames(CCLE_exp)[2:dim(CCLE_exp)[2]])

# replace gene to colnames
colnames(CCLE_exp)[2:dim(CCLE_exp)[2]] = genelist_exp

# get the paralog pair dependency 
paralog_sum_with_status_aneu = readRDS("paralog_sum_with_status_aneu.RDS")

# read the file all significant 37 significant paralogs
sig_dep_df = readxl::read_xlsx("dep_summary_aneuploidy_chr_loss.xlsx")
sig_dep_df %>% 
  dplyr::filter(status == "Paralog")-> sig_37_paralog_df

# get into paralog_pair_flt
# get the paralog_pairs
sig_37_paralog_df %>% 
  dplyr::mutate(para_cmb_1 = paste(dep_gene,paralog_gene, sep = "_")) %>% 
  dplyr::mutate(para_cmb_2 = paste(paralog_gene,dep_gene, sep = "_")) %>% 
  dplyr::select(para_cmb_1,para_cmb_2) -> sig_37_paralog_cmb

# get the query paralog_combination
sig_sig_para_list = c(sig_37_paralog_cmb$para_cmb_1, sig_37_paralog_cmb$para_cmb_2)

# check for the paralog family size and identical score 
paralog_pair_flt %>% 
  dplyr::mutate(condition = ifelse(para_name %in% sig_sig_para_list, "sig", "non_sig")) %>% 
  dplyr::mutate(condition = factor(condition, levels = c("sig", "non_sig"))) -> family_identity_score_df



############################################
#################Figure 2L #################
############################################
################################################
## gene expression difference in sig_paralogs ##
################################################
# get the copy number dataset
# take the cell line from CCLE cnv data
CCLE_gene_cnv = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /PortalOmicsCNGeneLog2.csv")

# cell and patient info
cell_line_id = read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/Cell_lines_annotations_20181226.txt")

# rename the gene name on CCLE_gene_CNV
# change the name of data
colnames(CCLE_gene_cnv)[1] = "DepMap_ID"
#dim(df_gene_dep)

# gene the gene name 
dim(CCLE_gene_cnv)
genelist_Dep = sub("\\.\\..*","",colnames(CCLE_gene_cnv)[2:dim(CCLE_gene_cnv)[2]])

# replace gene to colnames
colnames(CCLE_gene_cnv)[2:dim(CCLE_gene_cnv)[2]] = genelist_Dep
head(CCLE_gene_cnv[,1:6])

# get the Y loss cell line ids
results_chrY <- getBM(attributes = c("chromosome_name", "entrezgene_id", "hgnc_symbol"),
                      filters = "chromosome_name", values = "Y", mart = mart)

# alternative 
gene_info %>% 
  dplyr::filter(chromosome_name == "Y") ->results_chrY

Y_chr_genename = results_chrY$hgnc_symbol[results_chrY$hgnc_symbol != ""]

# get the table as Y chromosome cnv
Y_cnv = CCLE_gene_cnv[,c(1,which(colnames(CCLE_gene_cnv) %in% Y_chr_genename))]

# get the male cell line
cell_line_id %>% 
  dplyr::filter(Gender == "male") %>% 
  dplyr::select(depMapID) %>% 
  unlist(use.names = F) -> male_cellline

# get the male Y chromomr cnv df
Y_cnv %>% 
  dplyr::filter(DepMap_ID %in% male_cellline) -> Y_cnv_male

# Y chr loss sample
Y_loss_sample = Y_cnv_male[rowSums(is.na(Y_cnv_male))>25,]$DepMap_ID 
# Y chr normal sample
Y_normal = Y_cnv_male$DepMap_ID [!(Y_cnv_male$DepMap_ID %in% Y_loss_sample)]

# comparing the expression of 37 genes on corresponding chromosomes
# get the cell loc genes and its corrpsodning chromosome position
sig_37_paralog_df %>% 
  dplyr::select(aneuploid_loss_chr, paralog_gene) -> paralog_loss_chr_gene_df
# get the chromosome location and apply for the loop
chr_pos = paralog_loss_chr_gene_df$aneuploid_loss_chr
chr_pos_edt = c()

for (chr in chr_pos) {
  if(chr == "chrY"){
    chr_pos_edt = c(chr_pos_edt, chr)
  }else{
    chr_tmp = str_sub(chr, 4, nchar(chr))
    chr_edt = paste("X", str_split(chr_tmp, pattern = "_")[[1]][1],str_split(chr_tmp, pattern = "_")[[1]][2],sep = "" )
    chr_pos_edt = c(chr_pos_edt, chr_edt)
  }
}

# assign the new chromosome position 

paralog_loss_chr_gene_df$chr_pos = chr_pos_edt

# compare the relative expression level 
#Sig relative paralog gene expression 
sig_para_relative_exp = function(x){
  
  # input list 
  df_tpm = paralog_loss_chr_gene_df[x,]
  
  # input paralog gene 
  gene = df_tpm$paralog_gene
  
  # input chromsome pos
  chr_pos = df_tpm$chr_pos
  
  if(chr_pos == "chrY"){
    
    # get the expression data
    exp_flt_df = CCLE_exp[,c(1,which(colnames(CCLE_exp) == gene))]
    
    # loss exp 
    loss_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% Y_loss_sample, 2]
    
    # median loss exp
    median_loss_exp = median(loss_exp_list)
    
    # neutral exp
    neutral_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% Y_normal, 2]
    
    median_neutral_exp = median(neutral_exp_list)
    
    # normalised the expression
    relative_loss_exp = median_loss_exp -median_neutral_exp
    relative_neutral_exp = median_neutral_exp -median_neutral_exp
    
  }else{
    # get the corresponding aneuploidy condition 
    df_chr_info_df = df_aneu_score_flt[,c(1, which(colnames(df_aneu_score_flt) == chr_pos))]
    # chr loss aneuploid cells
    chr_loss_sample = df_chr_info_df[df_chr_info_df[,2] == -1, ]$X
    
    # chr_neutral aneuploid cells
    chr_neutral_sample = df_chr_info_df[df_chr_info_df[,2] == 0, ]$X
    
    # get the expression data
    exp_flt_df = CCLE_exp[,c(1,which(colnames(CCLE_exp) == gene))]
    
    # loss exp 
    loss_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% chr_loss_sample, 2]
    
    # median loss exp
    median_loss_exp = median(loss_exp_list)
    
    # neutral exp
    neutral_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% chr_neutral_sample, 2]
    
    median_neutral_exp = median(neutral_exp_list)
    
    # normalised the expression
    relative_loss_exp = median_loss_exp -median_neutral_exp
    relative_neutral_exp = median_neutral_exp -median_neutral_exp
    
  }
  
  df_return = data.frame(relative_loss_exp, relative_neutral_exp)
  colnames(df_return) = c("chr_loss", "chr_neutral")
  
  return(df_return)
} 

sog_37_paralog_rel_exp = lapply(c(1:dim(paralog_loss_chr_gene_df)[1]), sig_para_relative_exp)
sog_37_paralog_rel_exp_cmb = do.call(rbind,sog_37_paralog_rel_exp )
sog_37_paralog_rel_exp_cmb %>% 
  gather(chr_loss:chr_neutral, key = "condition", value = "value") -> sig_mdf_exp
sig_mdf_exp %>% 
  ggplot(aes(x = condition, y = value))+
  geom_boxplot(width = 0.25, 
               aes(color = condition, fill = condition),
               alpha = 0.5,
               outliers = F) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.15) +
  geom_rangeframe()+
  theme_tufte() +
  theme(legend.position = "none") +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  labs(x = "", y = "Rel. paralog gene expression (TPM)")+
  coord_cartesian(ylim = c(-1, 0.2)) -> P_sig_paralog_exp

# statistic analysis

#p = 2.73e-13
wilcox.test(sig_mdf_exp[sig_mdf_exp$condition == "chr_loss", ]$value, sig_mdf_exp[sig_mdf_exp$condition == "chr_neutral", ]$value)

############################################
#################Figure 2M #################
############################################

############################################
##   37 significant paralogs correlation  ##
############################################

sig_para_dep_exp_cor = function(x){
  
  # input list 
  df_tpm = sig_37_paralog_df[x,]
  
  # input paralog gene_dep
  gene_dep = df_tpm$dep_gene
  
  # input paralog_loss gene
  gene_loss = df_tpm$paralog_gene
  
  # get the expresssion of gene_loss
  exp_flt_df = CCLE_exp[,c(1,which(colnames(CCLE_exp) == gene_loss))]
  
  # get the dependency data
  dep_flt_list = df_gene_dep[,c(1, which((colnames(df_gene_dep) ==gene_dep )))]
  
  if(x == 31){
    exp_flt_df %>% 
      left_join(dep_flt_list, by = "DepMap_ID") %>% 
      dplyr::filter(.[[2]] > 0) %>% 
      dplyr::filter(DepMap_ID %in% aneuploid_cellline) -> df_flt_exp_dep
  }else{
    # combine the data
    exp_flt_df %>% 
      left_join(dep_flt_list, by = "DepMap_ID") %>% 
      dplyr::filter(.[[2]] > 1) %>% 
      dplyr::filter(DepMap_ID %in% aneuploid_cellline) -> df_flt_exp_dep
    
  }
  
  cor_score = cor.test(df_flt_exp_dep[,2],df_flt_exp_dep[,3])$estimate
  df_return =data.frame(gene_dep,gene_loss ,cor_score)
  return(df_return)
} 

# apply for all the sig paralogs 
cor_exp_dep_list = lapply(c(1:dim(sig_37_paralog_df)[1]), sig_para_dep_exp_cor)

cor_exp_dep_list_cmb = do.call(rbind,cor_exp_dep_list)

cor_exp_dep_list_cmb %>% 
  ggplot(aes(x = cor_score, y = ..scaled..)) +
  geom_density(fill = "#5271A8", color = "#5271A8", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray49") +
  #geom_rangeframe() +
  theme_tufte() +
  xlim(c(-1,1)) +
  labs(x = "Correlation", y = "Scaled density") 


############################################
#################Figure 2N #################
############################################
# identify the function of significants paralogs
sig_para_pairs = unique(c(sig_37_paralog_df$dep_gene, sig_37_paralog_df$paralog_gene))
library(clusterProfiler)
library(org.Hs.eg.db)
attributes <- c( "hgnc_symbol", "entrezgene_id")

# Fetch the gene information
gene_info_siglist <- getBM(attributes = attributes,
                           filters = "hgnc_symbol",
                           value = sig_para_pairs,
                           mart = mart,
                           useCache = T)

# GO analysis
ego <- enrichGO(gene          = sig_para_pairs,
                OrgDb         = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

View(ego[1:500,])

ego = as.data.frame(ego)
ego %>% 
  .[c(2,4,10,11,13,14, 18),] %>% 
  ggplot(aes(x = reorder(Description, -qvalue), y = -log10(qvalue))) +
  geom_point(aes(size = Count, color 
                 = qvalue), alpha = 0.75) +
  geom_rangeframe() +
  theme_tufte()+
  coord_flip()+
  scale_color_gradient2(low = "#E82035" ,mid = "#E5F3FD",  high = "#5371A9", midpoint = 5e-4)  +
  labs(x = "Description")
writexl::write_xlsx(ego, "sig_37_paralog_GO_analysis.xlsx")


############################################
#################Figure 2O #################
############################################
library(rbioapi)
########## question for the significant paralog pairs has more shared partner than non-significant one ############
# get the function to check for the shared protein and protein interaction partner in shared paralog partners 
CO_SHARED_INTERACTOR = function(paralog_pair){
  # get the paralog pair list
  
  co_shared_percentage_list = c()
  for (pair in paralog_pair) {
    # get the gene name
    gene_list = str_split(pair, pattern = "_")[[1]]
    
    # pair1
    gene_1 = gene_list[1]
    
    # pair2
    gene_2 = gene_list[2]
    
    # get into the interation genes
    
    # gene_1 interaction
    
    proteins_mapped_1 <- rba_string_map_ids(ids = gene_1,
                                            species = 9606)
    
    proteins_mapped_query = proteins_mapped_1$stringId
    if(!is.null(proteins_mapped_query)){
      int_partners_1 <- rba_string_interaction_partners(ids = proteins_mapped_query,
                                                        species = 9606,
                                                        required_score = 500,
                                                        limit = 100) 
      
      int_partners_1_sig = int_partners_1[int_partners_1$score >= 0.7,]
      
      gene_1_interaction = int_partners_1_sig$preferredName_B
      
    }
    
    
    # gene 2 interaction 
    proteins_mapped_2 <- rba_string_map_ids(ids = gene_2,
                                            species = 9606)
    
    proteins_mapped_query_2 = proteins_mapped_2$stringId
    if(!is.null(proteins_mapped_query_2)){
      int_partners_2 <- rba_string_interaction_partners(ids = proteins_mapped_query_2,
                                                        species = 9606,
                                                        required_score = 500,
                                                        limit = 100) 
      
      int_partners_2_sig = int_partners_2[int_partners_2$score >= 0.7,]
      
      gene_2_interaction = int_partners_2_sig$preferredName_B
      
    }
    
    # get the interaction gene
    intersection =intersect(gene_1_interaction, gene_2_interaction)
    
    # get the percentage
    
    shared_percentage = (length(intersection) / length(gene_1_interaction) + length(intersection) / length(gene_2_interaction)) / 2
    
    # add into the list
    co_shared_percentage_list = c(co_shared_percentage_list, shared_percentage)
    
  }
  
  return(co_shared_percentage_list)
}  
#### look the table of significant paralog pairs #####

sig_37_paralog_df %>% 
  dplyr::mutate(para_pair_gene = paste(dep_gene,paralog_gene, sep = "_")) -> sig_37_paralog_df_mdf

# get the info of co_shared interaction partner

# get th sig paralog shard interaction partner
sig_co_location = CO_SHARED_INTERACTOR(sig_37_paralog_df_mdf$para_pair_gene)
sig_co_location_df = data.frame(sig_co_location)
colnames(sig_co_location_df) = "value"
sig_co_location_df$condition = "sig"
# get the non-significant pairs
family_identity_score_df %>% 
  dplyr::filter(condition == "non_sig") -> non_sig_pairs

# get the non-sig interatcion partner
non_sig_co_location = CO_SHARED_INTERACTOR(non_sig_pairs$para_name)

# 

non_sig_co_location_df = data.frame(non_sig_co_location)
colnames(non_sig_co_location_df) = "value"
non_sig_co_location_df$condition = "non-sig"

cmb_df_co_interaction = rbind(sig_co_location_df, non_sig_co_location_df)

saveRDS(cmb_df_co_interaction, "sig_37_shared_interaction_rds")
library(ggthemes)
cmb_df_co_interaction %>% 
  dplyr::mutate(condition = factor(condition, 
                                   levels = c("sig", "non-sig"))) %>% 
  ggplot(aes(x = condition, y =value )) +
  geom_boxplot(
    aes(color = condition, fill = condition),
    width = 0.3,
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
  ) -> paired_interction


# test
# p-value 0.007333
wilcox.test(cmb_df_co_interaction[cmb_df_co_interaction$condition == "sig", ]$value, 
            cmb_df_co_interaction[cmb_df_co_interaction$condition == "non-sig", ]$value)


ggarrange(paired_interction, ncol = 2)

############################################
#################Figure 2P #################
############################################
colors = c("#EF8536","#3A76AF")

# identity score 
family_identity_score_df %>% 
  ggplot(aes(x = condition, y = mean_identical_score)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  #geom_jitter( aes(color = status),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte()+
  theme(legend.position = "None")+
  labs(x = "", y = "Identical score") +
  scale_color_manual(values = c(colors[1], "#7852A9"))+
  scale_fill_manual(values = c(colors[1], "#7852A9")) -> P_identity

# perform test
# p-value = 0.003676
t.test(family_identity_score_df[family_identity_score_df$condition == "sig", ]$mean_identical_score, 
       family_identity_score_df[family_identity_score_df$condition == "non_sig", ]$mean_identical_score)


############################################
#################Figure 2Q #################
############################################
# family_size
family_identity_score_df %>% 
  ggplot(aes(x = condition, y = family_size)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  #geom_jitter( aes(color = status),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte()+
  theme(legend.position = "None")+
  labs(x = "", y = "family_size ") +
  scale_color_manual(values = c(colors[1], "#7852A9"))+
  scale_fill_manual(values = c(colors[1], "#7852A9")) -> P_family

#p-value = 2.64e-06
t.test(family_identity_score_df[family_identity_score_df$condition == "sig", ]$family_size, 
       family_identity_score_df[family_identity_score_df$condition == "non_sig", ]$family_size)




############################################
#################Figure 2R #################
############################################
############################################################################################################
####       in 37 sig: # Are significant paralogs more highly-conserved than non-sig paralogs?           ####
############################################################################################################
library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

calc_gene_conservation <- function(
    genes,
    bigwig,
    region = c("exon","gene_body","promoter"),
    promoter_upstream = 2000,
    promoter_downstream = 200,
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    orgdb = org.Hs.eg.db
){
  region <- match.arg(region)
  
  # Map symbols -> ENTREZ
  map <- AnnotationDbi::select(orgdb, keys = genes, keytype = "SYMBOL", columns = "ENTREZID")
  map <- map[!is.na(map$ENTREZID), c("SYMBOL","ENTREZID")]
  map <- map[!duplicated(map$SYMBOL), ]
  rownames(map) <- map$ENTREZID
  if (!nrow(map)) stop("None of the provided symbols mapped to ENTREZ IDs.")
  
  # Build GRanges per gene
  if (region == "exon") {
    ranges_list <- GenomicFeatures::exonsBy(txdb, by = "gene")[map$ENTREZID]
  } else if (region == "gene_body") {
    genes_txdb <- GenomicFeatures::genes(txdb)
    ranges_list <- GenomicRanges::GRangesList(genes_txdb[map$ENTREZID])
  } else { # promoter
    tx_by_gene <- GenomicFeatures::transcriptsBy(txdb, "gene")[map$ENTREZID]
    pick <- lapply(tx_by_gene, function(gr) if (length(gr)==0) GRanges() else gr[which.max(width(gr))])
    pick <- GRangesList(pick)
    prom <- promoters(unlist(pick, use.names = FALSE),
                      upstream = promoter_upstream, downstream = promoter_downstream)
    split_idx <- rep(names(pick), lengths(pick))
    ranges_list <- split(prom, split_idx)
  }
  
  # Drop empties
  ok <- lengths(ranges_list) > 0
  ranges_list <- ranges_list[ok]
  if (length(ranges_list) == 0) stop("No genomic ranges found for the requested region type.")
  
  # Normalize seqlevel style to UCSC (chr1, chr2, ...)
  # Do this on the unlisted ranges, then split again.
  ul <- unlist(ranges_list, use.names = FALSE)
  seqlevelsStyle(ul) <- "UCSC"   # <-- SAFE: value is a character(1)
  ranges_list <- split(ul, rep(names(ranges_list), lengths(ranges_list)))
  
  # Helper: weighted stats
  wmean <- function(scores, widths) {
    if (!length(scores)) return(NA_real_)
    sum(scores * widths, na.rm = TRUE) / sum(widths[!is.na(scores)], na.rm = TRUE)
  }
  wmedian <- function(scores, widths) {
    if (!length(scores)) return(NA_real_)
    o <- order(scores)
    s <- scores[o]; w <- widths[o]
    cw <- cumsum(w)/sum(w)
    s[which.max(cw >= 0.5)]
  }
  
  # Prepare holders
  res_mean <- setNames(rep(NA_real_, length(ranges_list)), names(ranges_list))
  res_median <- res_mean
  res_nbases <- setNames(integer(length(ranges_list)), names(ranges_list))
  
  # Work chromosome-by-chromosome (based on the actual ranges)
  all_ranges <- unlist(ranges_list, use.names = FALSE)
  chrs <- unique(as.character(seqnames(all_ranges)))
  
  for (chr in chrs) {
    # pick all gene ranges on this chr
    chr_ranges <- lapply(ranges_list, function(gr) gr[as.character(seqnames(gr)) == chr])
    chr_ranges <- chr_ranges[lengths(chr_ranges) > 0]
    if (!length(chr_ranges)) next
    
    merged_chr <- reduce(unlist(GRangesList(chr_ranges), use.names = FALSE))
    bw_chr <- import(bigwig, which = merged_chr)
    if (length(bw_chr) == 0) next
    
    # compute per gene
    for (gene_id in names(chr_ranges)) {
      ol <- findOverlaps(chr_ranges[[gene_id]], bw_chr, ignore.strand = TRUE)
      if (!length(ol)) next
      q <- chr_ranges[[gene_id]][queryHits(ol)]
      s <- bw_chr[subjectHits(ol)]
      intersected <- pintersect(q, s, ignore.strand = TRUE)
      w <- width(intersected)
      scores <- mcols(s)$score
      res_mean[gene_id]   <- wmean(scores, w)
      res_median[gene_id] <- wmedian(scores, w)
      res_nbases[gene_id] <- sum(w)
    }
  }
  
  kept_dt <- data.frame(ENTREZID = names(ranges_list),
                        gene = map[names(ranges_list), "SYMBOL"], row.names = NULL)
  
  out <- data.frame(
    gene = kept_dt$gene,
    n_bases_scored = as.integer(res_nbases),
    mean_score = as.numeric(res_mean),
    median_score = as.numeric(res_median),
    stringsAsFactors = FALSE
  )
  out[!is.na(out$mean_score) & out$n_bases_scored > 0, ]
}


# Exon-level conservation using PhyloP
#non_sig_dep_genes_flt_2 = non_sig_dep_genes_flt[-which(non_sig_dep_genes_flt %in% "HMGB1P1")]
map_tbl <- AnnotationDbi::select(org.Hs.eg.db, 
                                 keys = c(sig_dep_genes, 
                                          non_sig_dep_genes_flt),
                                 keytype = "SYMBOL", 
                                 columns = "ENTREZID")
data_entriz_id = names(GenomicFeatures::exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene"))

# filter unmapped the genes
map_tbl %>% 
  dplyr::filter(ENTREZID %in% data_entriz_id) -> conser_gene_lists_flt

cons_all_paralog_exon <- calc_gene_conservation(
  genes = conser_gene_lists_flt$SYMBOL,
  bigwig = "~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /hg38.phyloP100way.bw",
  region = "exon"
)
head(cons_all_paralog_exon)
cons_all_paralog_exon %>% 
  dplyr::mutate(condition = ifelse(gene %in% sig_dep_genes, "sig", "non_sig" )) %>% 
  mutate(condition = factor(condition, levels = c("sig", "non_sig"))) -> conserv_df_plot
conserv_df_plot %>% 
  ggplot(aes(condition, median_score)) +
  
  geom_boxplot(
    aes(color = condition, fill = condition),
    width = 0.3,
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
  labs(x = "", y = "PhyloP score")+
  ylim(c(-1,2.5))  -> PhyloP_plot

#2.542e-05
wilcox.test(conserv_df_plot[conserv_df_plot$condition == "sig", ]$median_score, 
            conserv_df_plot[conserv_df_plot$condition == "non_sig", ]$median_score)

# plot the figures 
ggarrange(PhyloP_plot, ncol = 2)



############################################
#################Figure 2S #################
############################################
# Library
library(tidyverse)
library(ggthemes)
library(ggpubr)

# dataset
localisation_protein_df = read.delim("protein_localisation_dataset/subcellular_location.tsv")

# significant paralog
sig_paralog = readxl::read_xlsx("protein_localisation_dataset/sig_37_paralog.xlsx")

# non-significant paralog 
non_sig_paralog = readxl::read_xlsx("protein_localisation_dataset/non_sig_paralog.xlsx")

# first we have to know how many pair of can not find in the protein dataset in the protein atlas 
# Set rule: 
# either of paralog or neither of paralog gene in the dataset , we will exclude the analysis afterwards

# function 
# index here as input data
OVERLAPPED_FREQ = function(data, index_length){
  
  # generate vector
  frequency =  c()
  # get the input data
  for (i in c(1:index_length)) {
    # get the data by row
    df_tpm = data[i,]
    
    # get the paralog_1
    para_1 = unlist(df_tpm[,1],use.names = F)
    
    # get the paralog_2
    para_2 = unlist(df_tpm[,2],use.names = F)
    
    # look in to the location data sheet and get the location of protein (approved or supported)
    # both conditions are included into our analysis 
    
    # location of paralog1
    localisation_protein_df %>% 
      dplyr::filter(Gene.name == para_1) -> loc_flt_para_1
    
    # get the main positionand additional position 
    position_para_1 = c(loc_flt_para_1$Main.location,loc_flt_para_1$Additional.location)
    position_para_1_splt = unlist(strsplit(position_para_1, ";"))
    
    # location of paralog2
    localisation_protein_df %>% 
      dplyr::filter(Gene.name == para_2) -> loc_flt_para_2
    
    # get the main positionand additional position 
    position_para_2 = c(loc_flt_para_2$Main.location,loc_flt_para_2$Additional.location)
    position_para_2_splt = unlist(strsplit(position_para_2, ";"))
    
    # check for the intersection 
    overlapped_region = intersect(position_para_1_splt, position_para_2_splt)
    
    # calculate the mean frequency
    mean_freq = (length(overlapped_region)/length(position_para_1_splt) + 
                   length(overlapped_region)/length(position_para_2_splt)) / 2
    
    frequency = c(frequency, mean_freq)  
  }
  
  return(frequency)
}

# first check for the significant paralog pairs
sig_count = 0
sig_position = c()
for (i in c(1: dim(sig_paralog)[1])) {
  
  # get the df
  df_tmp = sig_paralog[i,]
  # get the paralog gene
  para_1 = unlist(df_tmp[,1],use.names = F)
  
  # get the paralog partner 
  para_2 = unlist(df_tmp[,7],use.names = F)
  
  # check if both of gene in the database of protein localation dataset
  # check for the parlaog 1
  sum_1 = sum(para_1 %in% localisation_protein_df$Gene.name)
  
  # check for the paralog 2
  sum_2 = sum(para_2 %in% localisation_protein_df$Gene.name)
  
  # if sum1 + sum2 = 2 means both sample are presenting in the dataset 
  if(sum_1 + sum_2 == 2){
    sig_count = sig_count + 1
    
    # add into the position list 
    sig_position = c(sig_position, i)
  }
  
  
}

# check for the non-significant 
non_sig_count = 0
non_sig_position = c()
for (i in c(1: dim(non_sig_paralog)[1])) {
  
  # get the df
  df_tmp = non_sig_paralog[i,]
  # get the paralog gene
  para_1 = unlist(df_tmp[,3],use.names = F)
  
  # get the paralog partner 
  para_2 = unlist(df_tmp[,4],use.names = F)
  
  # check if both of gene in the database of protein localation dataset
  # check for the parlaog 1
  sum_1 = sum(para_1 %in% localisation_protein_df$Gene.name)
  
  # check for the paralog 2
  sum_2 = sum(para_2 %in% localisation_protein_df$Gene.name)
  
  # if sum1 + sum2 = 2 means both sample are presenting in the dataset 
  if(sum_1 + sum_2 == 2){
    non_sig_count = non_sig_count + 1
    
    # add into the position list 
    non_sig_position = c(non_sig_position, i)
  }
  
  
}

# plot the presenting of sig and non-sig paralog pair presenting in the dataset 
# 25/37
# 13062/28103
summary_df = tibble::tibble(freq = c(sig_count/dim(sig_paralog)[1], non_sig_count/dim(non_sig_paralog)[1]),
                            annotation = c("significant paralog", "non-significant paralog"))

# plot the statistic 
summary_df %>% 
  ggplot(aes(x = annotation, y = freq)) +
  geom_col(aes(color = annotation, 
               fill = annotation), width = 0.5, alpha = 0.7) +
  theme_tufte() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# check for the position that is significant 
sig_flt_df = sig_paralog[sig_position,c(1,7)]


# to check for the overlapping frequency in the significant pairs
freq_sig = OVERLAPPED_FREQ(sig_flt_df, dim(sig_flt_df)[1]) 
freq_df_sig = data.frame(freq_sig)
freq_df_sig$annotation = "sig"
colnames(freq_df_sig)[1] = "freq"
# to check for the non-sig freq
nonsig_flt_df = non_sig_paralog[non_sig_position,c(3,4)]

# to check for the overlapping frequency in the non_significant pairs
freq_non_sig = OVERLAPPED_FREQ(nonsig_flt_df, dim(nonsig_flt_df)[1]) 

freq_df_nonsig = data.frame(freq_non_sig)
freq_df_nonsig$annotation = "non-sig"
colnames(freq_df_nonsig)[1] = "freq"

# combine the data
freq_df_cmb = rbind(freq_df_sig, freq_df_nonsig)

colors = c("#EF8536","#3A76AF")
# plot the figure
freq_df_cmb %>% 
  mutate(annotation = factor(annotation, levels = c("sig", "non-sig"))) %>% 
  ggplot(aes(x = annotation, y = freq)) +
  geom_boxplot(width = 0.25, 
               aes(color = annotation, fill = annotation),
               alpha = 0.5,
               outliers = F) +
  #geom_jitter( aes(color = annotation),position=position_jitter(0.12),size = 1,alpha = 0.15) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[1], "#7852A9"))+
  scale_fill_manual(values = c(colors[1], "#7852A9"))

####### t.test ##########
# p-value = 0.01451
t.test(freq_df_cmb[freq_df_cmb$annotation == "sig",]$freq, 
       freq_df_cmb[freq_df_cmb$annotation == "non-sig",]$freq)

