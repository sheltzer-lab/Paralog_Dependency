###################################################
########### Figure S2 code ########################
###################################################
# Author: Yi
# Date: Oct_2026

# project :calculate the VAF from the sequecing data
# Date: April.16.2025
# Author: Yi

# load library
library(tidyverse)
library(ggthemes)
library(ggpubr)

# load data
# list file in VAF folder
File_list = list.files("VAF_folder",full.names = T)
sample_name = c()
gene_list = c()
for (i in File_list) {
  
  # read each line of the input
  file_name = i
  info_tpm = str_split(file_name, pattern = "/")[[1]][2]
  
  # split the string and get the sample name and gene name 
  gene_name = str_split(info_tpm,pattern = "_")[[1]][[2]]
  sample_name_tpm = paste(str_split(info_tpm,pattern = "_")[[1]][[1]], 
                          str_split(info_tpm,pattern = "_")[[1]][[2]] ,
                          str_split(info_tpm,pattern = "_")[[1]][[3]],sep = "_"
  )
  
  sample_name = c(sample_name,sample_name_tpm)
  gene_list = c(gene_list, gene_name)
}


# make function and read all the VAF files 
READ_VAF = function(x){
  file_name = File_list[x]
  Sample_ID = sample_name[x]
  gene = gene_list[x]
  lines <- readLines(file_name)
  split_lines <- strsplit(lines, "\t")
  parsed <- lapply(split_lines, function(x) {
    chrom <- x[1]
    pos <- as.integer(x[2])
    ref <- x[3]
    depth <- as.integer(x[4])
    
    # Process base fields (can be A, C, G, T, DEL, etc.)
    base_data <- x[5:length(x)]
    base_df <- do.call(rbind, lapply(base_data, function(b) {
      parts <- unlist(strsplit(b, ":"))
      data.frame(
        base = parts[1],
        count = as.integer(parts[2]),
        avg_mq = as.numeric(parts[3]),
        avg_bq = as.numeric(parts[4]),
        stringsAsFactors = FALSE
      )
    }))
    
    base_df$chrom <- chrom
    base_df$pos <- pos
    base_df$ref <- ref
    base_df$depth <- depth
    base_df$sample <- Sample_ID
    base_df$gene <- gene
    
    return(base_df)
  })
  
  result_df <- bind_rows(parsed)
  return(result_df)
}

VAF_list = lapply(c(1:length(File_list)), READ_VAF)
VAF_list_cmb = do.call(rbind, VAF_list)

#write_rds(VAF_list_cmb, "VAF_list_cmb.RDS")
VAF_list_cmb = readRDS("VAF_list_cmb.RDS")
# read_sample sheet with position and gene and guide order
library(readxl)
sgRNA_datasheet = read_xlsx("vaf_sgRNA_query_position.xlsx")

# calculate the VAF in the files 
VAF_calculation = function(x){
  # get the sample name
  sample_id = sample_name[x]
  gene = gene_list[x]
  # sgRNA guide 
  guide = paste(str_split(sample_id, pattern = "_")[[1]][2], str_split(sample_id, pattern = "_")[[1]][3], sep = "_")
  
  # get the position 
  sgRNA_datasheet %>% 
    filter(target_order == guide) ->df_pos
  
  # get the chr nad pos
  chr = df_pos$chr
  pos = df_pos$position
  
  # check the region [-1, 0, 1]
  
  check_pos = c(-10:10) # base pair interval
  mutant_count_list = c()
  total_read = c()
  if(!(gene %in% c("sgAAVS1", "sgNT","sgRosa26"))){
    for (i in check_pos) {
      pos_tpm = pos + i
      # get the sample after filter and get the query data
      VAF_list_cmb %>% 
        dplyr::filter(sample == sample_id) %>% 
        dplyr::filter(chrom == chr,
                      pos == pos_tpm) -> df_tpm
      
      
      #get the ref target
      ref = unique(df_tpm$ref)
      
      # get the total count
      total_count = unique(df_tpm$depth)
      
      # get the mutation read
      df_tpm %>% 
        dplyr::filter(base != ref) -> df_tpm_flt
      
      mutant_count = sum(df_tpm_flt$count)
      mutant_count_list = c(mutant_count_list ,mutant_count )
      total_read = c(total_read, total_count)
    }
    
    VAF_freq = sum(mutant_count_list) / mean(total_read)
    df_return = data.frame(gene,sample_id, guide, VAF_freq)
  }else{
    VAF_freq = NA
    df_return = data.frame(gene,sample_id, guide, VAF_freq)
  }
  
  return(df_return)
  
}

VAF_cal_guide = lapply(c(1:length(sample_name)), VAF_calculation)
VAF_cal_guide_cmb = do.call(rbind, VAF_cal_guide)

VAF_cal_guide_cmb %>% 
  dplyr::filter(!is.na(VAF_freq)) -> target_sample

# Cell line A
target_sample %>% 
  dplyr::mutate(sample_id2 = sample_id) %>% 
  separate(col = sample_id2, sep = "_", c("cell_line","gene2","index")) %>%
  dplyr::filter(cell_line =="A") %>% 
  ggplot(aes(x = guide, y= VAF_freq))+
  geom_col(aes(color = gene,
               fill = gene),
           width = 0.5,
           alpha = 0.8)+
  geom_hline(yintercept = 0.05,linetype = "dashed",color="red")+
  geom_rangeframe()+
  theme_tufte() +
  theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  scale_color_manual(values = c("#5271A8","#6895BD",  "#7852A9","#dddddd")) +
  scale_fill_manual(values = c("#5271A8","#6895BD",  "#7852A9","#dddddd")) +
  #facet_wrap(~cell_line,nrow = 1, ncol = 3)+
  labs(x = "guides", y = "VAF")+
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 0.7)) + 
  ggtitle(label = "A") -> F_A

# Cell line H
target_sample %>% 
  dplyr::mutate(sample_id2 = sample_id) %>% 
  separate(col = sample_id2, sep = "_", c("cell_line","gene2","index")) %>%
  dplyr::filter(cell_line =="H") %>% 
  ggplot(aes(x = guide, y= VAF_freq))+
  geom_col(aes(color = gene,
               fill = gene),
           width = 0.5,
           alpha = 0.8)+
  geom_hline(yintercept = 0.05,linetype = "dashed",color="red")+
  geom_rangeframe()+
  theme_tufte() +
  theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  scale_color_manual(values = c("#5271A8","#6895BD",  "#7852A9","#dddddd")) +
  scale_fill_manual(values = c("#5271A8","#6895BD",  "#7852A9","#dddddd")) +
  #facet_wrap(~cell_line,nrow = 1, ncol = 3)+
  labs(x = "guides", y = "VAF") +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 0.7)) + 
  ggtitle(label = "H") -> F_H

# Cell line M
target_sample %>% 
  dplyr::mutate(sample_id2 = sample_id) %>% 
  separate(col = sample_id2, sep = "_", c("cell_line","gene2","index")) %>%
  dplyr::filter(cell_line =="M") %>%
  dplyr::filter(!(guide %in% c("PPP2CA_1", "PPP2CA_2", "PPP2CA_4"))) %>% 
  ggplot(aes(x = guide, y= VAF_freq))+
  geom_col(aes(color = gene,
               fill = gene),
           width = 0.5,
           alpha = 0.8)+
  geom_hline(yintercept = 0.05,linetype = "dashed",color="red")+
  geom_rangeframe()+
  theme_tufte() +
  theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  scale_color_manual(values = c("#5271A8","#6895BD",  "#7852A9","#dddddd")) +
  scale_fill_manual(values = c("#5271A8","#6895BD",  "#7852A9","#dddddd")) +
  #facet_wrap(~cell_line,nrow = 1, ncol = 3)+
  labs(x = "guides", y = "VAF") +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 0.7)) + 
  ggtitle(label = "M") -> F_M

ggarrange(F_A, F_H,F_M,
          nrow = 1, ncol = 3, align = "hv")

# check the control guides with editing frequency
# only focus on AAVS1_1, AAVS1_2, sgNT and sgRosa26
#
control_sample = c("A_sgAAVS1_1", "A_sgAAVS1_2", "A_sgNT_vaf", "A_sgRosa26_vaf",
                   "H_sgAAVS1_1", "H_sgAAVS1_2", "H_sgNT_vaf", "H_sgRosa26_vaf",
                   "M_sgAAVS1_1", "M_sgAAVS1_2", "M_sgNT_vaf", "M_sgRosa26_vaf")

VAF_calculation_control_guide = function(x){
  # get the sample name
  sample_id = control_sample[x]
  
  guides_list = c()
  VAF_list = c()
  
  # sgRNA guide 
  for (i in c(1:dim(sgRNA_datasheet)[1])) {
    guides = sgRNA_datasheet[i,]$target_order
    
    # # get the chr nad pos
    chr = sgRNA_datasheet[i,]$chr
    pos = sgRNA_datasheet[i,]$position
    
    check_pos = c(-10:10)
    
    
    mutant_count_list = c()
    total_read = c()
    
    
    for (j in check_pos) {
      pos_tpm = pos + j
      # get the sample after filter and get the query data
      VAF_list_cmb %>% 
        dplyr::filter(sample == sample_id) %>% 
        dplyr::filter(chrom == chr,
                      pos == pos_tpm) -> df_tpm
      
      
      #get the ref target
      ref = unique(df_tpm$ref)
      
      # get the total count
      total_count = unique(df_tpm$depth)
      
      # get the mutation read
      df_tpm %>% 
        dplyr::filter(base != ref) -> df_tpm_flt
      
      mutant_count = sum(df_tpm_flt$count)
      mutant_count_list = c(mutant_count_list ,mutant_count )
      total_read = c(total_read, total_count)
    }
    
    VAF_freq = sum(mutant_count_list) / mean(total_read)
    
    guides_list = c(guides_list,guides)
    VAF_list = c(VAF_list,VAF_freq)
    
  }
  
  return_df = data.frame(guides_list, VAF_list,sample_id)
  return(return_df)
}

# get the frequency of 
control_mutation_vaf = lapply(c(1: length(control_sample)), VAF_calculation_control_guide)
control_mutation_vaf_cmb = do.call(rbind, control_mutation_vaf)
# get the sgRNA order
guides_order = target_sample[1:15,]$guide

control_mutation_vaf_cmb %>% 
  dplyr::mutate(guides_list = factor(guides_list, levels = guides_order)) %>% 
  filter(!is.na(guides_list)) %>% 
  dplyr::mutate(guides_list2 = guides_list) %>% 
  separate(col = guides_list2, sep = "_", c("gene","slash")) %>% 
  ggplot(aes(x = guides_list, y= VAF_list)) +
  geom_col(aes(color = gene,
               fill = gene),width = 0.5)+
  geom_hline(yintercept = 0.05,linetype = "dashed",color = "red")+
  geom_rangeframe() +
  theme_tufte() +
  theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  scale_color_manual(values = c("#5271A8","#6895BD",  "#7852A9","#dddddd")) +
  scale_fill_manual(values = c("#5271A8","#6895BD",  "#7852A9","#dddddd")) +
  facet_wrap(~sample_id) +
  labs(x = "guides", y = "VAF")

# check for its corresponding paralog off-target frequency
# gene: PPP2CB, PPP2R1B, DHFR, BTRC
# the gene position :
# PPP2CB: chr8 30785616-30812818
# PPP2R1B: chr11 111688000-111766389
# DHFR: chr5 80626226 - 80654983
# BTRC: chr10 101354048- 101557313

# let us look into each paralog genes and check for the general frequency of mutation 
paralog_query_list = tibble("chr" = c("chr8", "chr11", "chr5", "chr10"),
                            "start" = c(30785616, 111688000, 80626226, 101354048),
                            "end" = c(30812818, 111766389, 80654983, 101557313),
                            "gene" = c("PPP2CB", "PPP2R1B", "DHFR","BTRC"))
paralog_query_list = as.data.frame(paralog_query_list)

offtarget_paralog_freq = function(x){
  # get sample name
  sample_id = sample_name[x]
  
  # get the data for 
  VAF_list_cmb %>% 
    dplyr::filter(sample == sample_id) -> df_tpm
  
  # generate the list
  gene_list = c()
  VAF_list = c()
  
  # get the list of gene 
  for (i in c(1:dim(paralog_query_list)[1])) {
    # gene name
    gene_name = paralog_query_list[i,]$gene
    
    # start site
    start_site = paralog_query_list[i,]$start
    
    # end site
    end_site = paralog_query_list[i,]$end
    # get the range of gene locus
    df_tpm %>% 
      dplyr::filter(pos >= start_site & pos <= end_site ) -> df_tpm_flt
    
    mutation_list = c()
    total_read = c()
    # get the position size
    pos_list = unique(df_tpm_flt$pos)
    for (i in pos_list) {
      
      df_tpm_flt %>% 
        dplyr::filter(pos == i) -> df_tpm_flt_site
      
      # get the ref 
      ref = df_tpm_flt_site$ref
      
      df_tpm_flt_site %>% 
        filter(base != ref) -> df_tpm_flt_ref
      
      # get the sum of mutant reads
      sum_mutation = sum(df_tpm_flt_ref$count)
      mutation_list = c(mutation_list, sum_mutation)
      # get the total reads 
      total = df_tpm_flt_ref$depth
      total_read = c(total_read,total)
    }
    
    # calculate the frequency of VAF 
    VAF = sum(mutation_list) / sum(total_read)
    
    # add to the list 
    gene_list = c(gene_list, gene_name)
    VAF_list = c(VAF_list, VAF)
  }
  
  return_df = data.frame(sample_id, gene_list, VAF_list)
  return(return_df)
}

paralog_offtarget_cal = lapply(c(1:length(sample_name)), offtarget_paralog_freq)

paralog_offtarget_cal_cmb = do.call(rbind, paralog_offtarget_cal)
#write_rds(paralog_offtarget_cal_cmb,"paralog_offtarget_cal_cmb.rds")
paralog_offtarget_cal_cmb = readRDS("paralog_offtarget_cal_cmb.rds")

# get the cell line name
cell_line_list = c()
# get the ontarget gene name
ontarget_gene = c()
#guide list
guides_list = c()
for (i in paralog_offtarget_cal_cmb$sample_id) {
  vector_tpm = i
  cellname = str_split(vector_tpm,pattern = "_")[[1]][1]
  ontarget = str_split(vector_tpm,pattern = "_")[[1]][2]
  guides = paste(str_split(vector_tpm,pattern = "_")[[1]][2], str_split(vector_tpm,pattern = "_")[[1]][3],sep = "_")
  # add to the list
  cell_line_list = c(cell_line_list,cellname )
  ontarget_gene = c(ontarget_gene,ontarget )
  guides_list = c(guides_list, guides)
}

# get the ontarget gene name
paralog_offtarget_cal_cmb$cell = cell_line_list
paralog_offtarget_cal_cmb$ontarget = ontarget_gene
paralog_offtarget_cal_cmb$guides = guides_list
# PPP2CB
paralog_offtarget_cal_cmb %>% 
  dplyr::filter(gene_list == "PPP2CB") %>% 
  dplyr::filter(ontarget %in% c("PPP2CA","sgAAVS1","sgNT","sgRosa26")) %>% 
  ggplot(aes(x = guides ,y = VAF_list ))+
  geom_col(aes(color = ontarget,
               fill = ontarget),
           width = 0.6, alpha = 0.75) +
  geom_rangeframe()+
  theme_tufte()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~cell, nrow = 3, ncol = 1) +
  scale_color_manual(values = c("#6895BD","#FFE5E1",  "#FFC1B9","#FF4F37")) +
  scale_fill_manual(values = c("#6895BD","#FFE5E1",  "#FFC1B9","#FF4F37")) +
  ggtitle("PPP2CB") +
  labs(x = "Guides", y = "VAF") -> p1

#PPP2R1B
paralog_offtarget_cal_cmb %>% 
  dplyr::filter(gene_list == "PPP2R1B") %>% 
  dplyr::filter(ontarget %in% c("PPP2R1A","sgAAVS1","sgNT","sgRosa26")) %>% 
  ggplot(aes(x = guides ,y = VAF_list ))+
  geom_col(aes(color = ontarget,
               fill = ontarget),
           width = 0.6, alpha = 0.75) +
  geom_rangeframe()+
  theme_tufte()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~cell, nrow = 3, ncol = 1) + 
  scale_color_manual(values = c("#7852A9", "#FFE5E1",  "#FFC1B9","#FF4F37")) +
  scale_fill_manual(values = c("#7852A9",  "#FFE5E1",  "#FFC1B9","#FF4F37")) +
  ggtitle("PPP2R1B") +
  labs(x = "Guides", y = "VAF") -> p2

#
paralog_offtarget_cal_cmb %>% 
  dplyr::filter(gene_list == "DHFR") %>% 
  dplyr::filter(ontarget %in% c("TYMS","sgAAVS1","sgNT","sgRosa26")) %>% 
  dplyr::mutate(ontarget = factor(ontarget, levels = c("TYMS", "sgAAVS1","sgNT","sgRosa26"))) %>% 
  ggplot(aes(x = factor(guides,
                        levels =c("TYMS_1","TYMS_2","TYMS_3","TYMS_4",
                                  "sgAAVS1_1","sgAAVS1_2","sgNT_vaf","sgRosa26_vaf")) ,
             y = VAF_list ))+
  geom_col(aes(color = ontarget,
               fill = ontarget),
           width = 0.6, alpha = 0.75) +
  geom_rangeframe()+
  theme_tufte()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~cell, nrow = 3, ncol = 1) +
  scale_color_manual(values = c( "#dddddd","#FFE5E1",  "#FFC1B9","#FF4F37")) +
  scale_fill_manual(values = c( "#dddddd","#FFE5E1",  "#FFC1B9","#FF4F37")) +
  ggtitle("DHFR") +
  labs(x = "Guides", y = "VAF")  -> p3

#BTRC
paralog_offtarget_cal_cmb %>% 
  dplyr::filter(gene_list == "BTRC") %>% 
  dplyr::filter(ontarget %in% c("FBXW11","sgAAVS1","sgNT","sgRosa26")) %>% 
  ggplot(aes(x = guides ,y = VAF_list ))+
  geom_col(aes(color = ontarget,
               fill = ontarget),
           width = 0.6, alpha = 0.75) +
  geom_rangeframe()+
  theme_tufte()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~cell, nrow = 3, ncol = 1) +
  scale_color_manual(values = c("#5271A8","#FFE5E1",  "#FFC1B9","#FF4F37")) +
  scale_fill_manual(values = c("#5271A8","#FFE5E1",  "#FFC1B9","#FF4F37")) +
  ggtitle("BTRC") +
  labs(x = "Guides", y = "VAF") -> p4

ggarrange(p4,p1,p2,p3,
          nrow = 1,ncol = 4)

# calculate the p value 
cell_types = unique(paralog_offtarget_cal_cmb$cell)
gene_list = unique(paralog_offtarget_cal_cmb$ontarget)
# statistic 
cell_lists = c()
ontarget_lists = c()
offtarget_list = c()
pvalue_lists = c()
for(i in cell_types){
  for (j in gene_list) {
    cell_name = i
    gene_name = j
    
    # select
    if(gene_name == "FBXW11"){
      paralog_offtarget_cal_cmb %>% 
        dplyr::filter(gene_list == "BTRC") %>% 
        dplyr::filter(cell == cell_name,
                      ontarget %in% c(gene_name, "sgAAVS1","sgNT","sgRosa26") ) ->df_tpm
      
      # p_value
      p_value = t.test(df_tpm[df_tpm$ontarget == gene_name, ]$VAF_list,df_tpm[df_tpm$ontarget %in% c("sgAAVS1","sgNT","sgRosa26"), ]$VAF_list )$p.value
      cell_lists= c(cell_lists, cell_name)
      ontarget_lists = c(ontarget_lists, gene_name)
      offtarget_list = c(offtarget_list , "BTRC")
      pvalue_lists = c(pvalue_lists, p_value)
      
    }
    
    if(gene_name == "PPP2CA"){
      paralog_offtarget_cal_cmb %>% 
        dplyr::filter(gene_list == "PPP2CB") %>% 
        dplyr::filter(cell == cell_name,
                      ontarget %in% c(gene_name, "sgAAVS1","sgNT","sgRosa26") )->df_tpm
      
      # p_value
      p_value = t.test(df_tpm[df_tpm$ontarget == gene_name, ]$VAF_list,df_tpm[df_tpm$ontarget %in% c("sgAAVS1","sgNT","sgRosa26"), ]$VAF_list )$p.value
      cell_lists= c(cell_lists, cell_name)
      ontarget_lists = c(ontarget_lists, gene_name)
      offtarget_list = c(offtarget_list , "PPP2CB")
      pvalue_lists = c(pvalue_lists, p_value)
    }
    
    if(gene_name == "PPP2R1A"){
      paralog_offtarget_cal_cmb %>% 
        dplyr::filter(gene_list == "PPP2R1B") %>% 
        dplyr::filter(cell == cell_name,
                      ontarget %in% c(gene_name, "sgAAVS1","sgNT","sgRosa26") )->df_tpm
      
      # p_value
      p_value = t.test(df_tpm[df_tpm$ontarget == gene_name, ]$VAF_list,df_tpm[df_tpm$ontarget %in% c("sgAAVS1","sgNT","sgRosa26"), ]$VAF_list )$p.value
      cell_lists= c(cell_lists, cell_name)
      ontarget_lists = c(ontarget_lists, gene_name)
      offtarget_list = c(offtarget_list , "PPP2R1B")
      pvalue_lists = c(pvalue_lists, p_value)
    }
    
    if(gene_name == "TYMS"){
      paralog_offtarget_cal_cmb %>% 
        dplyr::filter(gene_list == "DHFR") %>% 
        dplyr::filter(cell == cell_name,
                      ontarget %in% c(gene_name, "sgAAVS1","sgNT","sgRosa26") )->df_tpm
      
      
      # p_value
      p_value = t.test(df_tpm[df_tpm$ontarget == gene_name, ]$VAF_list,df_tpm[df_tpm$ontarget %in% c("sgAAVS1","sgNT","sgRosa26"), ]$VAF_list )$p.value
      cell_lists= c(cell_lists, cell_name)
      ontarget_lists = c(ontarget_lists, gene_name)
      offtarget_list = c(offtarget_list , "DHFR")
      pvalue_lists = c(pvalue_lists, p_value)
    }
    
  }
}

stat_summary = data.frame(cell_lists,ontarget_lists ,offtarget_list,pvalue_lists)
# check for the vcf file and call and calculate the TMB
# get the bed file

library(data.table)
bed <- fread("S07604514_Regions.bed", header = FALSE)
colnames(bed) <- c("chrom", "start", "end")
bed %>% 
  dplyr::mutate(diff = (end- start) +1) -> bed_edt
# the coverage is 38.4493 Mb
total_region = sum(bed_edt$diff)/1000000


#run VCF files
vcf_files = list.files("/Volumes/My Passport for Mac/Mutect2_annotate_vcf",pattern = "\\.vcf$",full.names = T)[1:57]
library(VariantAnnotation)
sample_names = c()
cell_lines = c()
gene_lists = c()
guides =c()
for (i in vcf_files) {
  file_tpm = str_split(i, pattern = "/")[[1]][5]
  sample_name = paste(str_split(file_tpm, pattern = "_")[[1]][1],
                      str_split(file_tpm, pattern = "_")[[1]][2],
                      str_split(file_tpm, pattern = "_")[[1]][3],sep = "_")
  cell_line = str_split(file_tpm, pattern = "_")[[1]][1]
  gene = str_split(file_tpm, pattern = "_")[[1]][2]
  guide =  paste(
    str_split(file_tpm, pattern = "_")[[1]][2],
    str_split(file_tpm, pattern = "_")[[1]][3],sep = "_")
  # add to the vector
  sample_names = c(sample_names, sample_name)
  cell_lines = c(cell_lines, cell_line)
  gene_lists = c(gene_lists, gene)
  guides = c(guides, guide)
}

TMB_cal = function(x){
  # read file
  df_tpm = readVcf(vcf_files[x],genome = "hg38")
  
  sample = sample_names[x]
  cell_line = cell_lines[x]
  gene = gene_lists[x]
  guide = guides[x]
  # Keep only variants that passed all filters
  vcf_pass <- df_tpm[fixed(df_tpm)$FILTER == "PASS", ]
  
  # only keep SNV
  is_snv <- nchar(as.character(ref(vcf_pass))) == 1 & sapply(alt(vcf_pass), function(x) all(nchar(as.character(x)) == 1))
  vcf_snv <- vcf_pass[is_snv, ]
  
  
  tmb_count <- nrow(vcf_snv)
  tmb <- tmb_count / total_region  # result: mutations per Mb
  return_df = data.frame(sample, cell_line, gene, guide, tmb)
  
  return(return_df)
}

TMB_sum = lapply(c(1:length(vcf_files)), TMB_cal)
TMB_sum_cmb = do.call(rbind, TMB_sum)

# plot
TMB_sum_cmb %>% 
  dplyr::mutate(gene = factor(gene, levels = c("FBXW11","PPP2CA","PPP2R1A","TYMS", "sgAAVS1","sgNT","sgRosa26"))) %>%
  ggplot(aes(x = factor(guide, levels = c("FBXW11_1","FBXW11_3","FBXW11_4",
                                          "PPP2CA_1","PPP2CA_2","PPP2CA_3","PPP2CA_4",
                                          "PPP2R1A_1","PPP2R1A_2","PPP2R1A_3","PPP2R1A_4",
                                          "TYMS_1","TYMS_2","TYMS_3","TYMS_4",
                                          "sgAAVS1_1","sgAAVS1_2","sgNT_annotated","sgRosa26_annotated")),
             y =tmb))+
  geom_col(aes(color = gene,
               fill = gene),
           width = 0.6, alpha = 0.75)+
  geom_rangeframe()+
  theme_tufte()+
  scale_color_manual(values = c("#5271A8","#6895BD",  "#7852A9","#dddddd","#FFE5E1",  "#FFC1B9","#FF4F37")) +
  scale_fill_manual(values = c("#5271A8","#6895BD",  "#7852A9","#dddddd","#FFE5E1",  "#FFC1B9","#FF4F37")) +
  facet_wrap(~cell_line,nrow = 1,ncol = 3,scales = "free_y")+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  labs(y = "TMB")
