################################
########### Figure 3 ###########
################################
# Author: Yi
# Date: Oct_2026

# Project: all protein coding gene TSS in chromosome loss and neutral cell (all chromosome interagtion)
# Author: Yi Bei
# Date: Oct.02.2026

# library
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)


############################################
#################Figure 3B #################
############################################
#data
sample_atac = read.delim("../TCGA_identifier_mapping.txt")
TACG_df = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Aneuploidy patient datasets/mmc2.xlsx", skip = 1)


# function 
# ----------------------------
# function to read one matrix
# ----------------------------
read_profile_tss <- function(file) {
  
  x <- readLines(file)
  
  # line 1: group names and group sizes
  h1 <- sub("^#", "", x[1])
  ginfo <- strsplit(h1, "\t")[[1]]
  
  gnames <- sub(":(.*)$", "", ginfo)
  gsizes <- as.numeric(sub("^.*:", "", ginfo))
  
  # line 2: parameters
  h2 <- sub("^#", "", x[2])
  
  upstream   <- as.numeric(sub(".*upstream:(\\d+).*", "\\1", h2))
  downstream <- as.numeric(sub(".*downstream:(\\d+).*", "\\1", h2))
  bin_size   <- as.numeric(sub(".*bin size:(\\d+).*", "\\1", h2))
  
  # numeric matrix starts after first 3 lines
  num_lines <- x[-c(1, 2, 3)]
  
  mat <- do.call(
    rbind,
    lapply(num_lines, function(z) as.numeric(strsplit(z, "\t")[[1]]))
  )
  
  mat <- as.matrix(mat)
  
  group <- rep(gnames, gsizes)
  
  out <- do.call(
    rbind,
    lapply(gnames, function(g) {
      idx <- which(group == g)
      
      data.frame(
        group = g,
        bin = seq_len(ncol(mat)),
        signal = colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
      )
    })
  )
  
  out$sample <- basename(file)
  out$upstream <- upstream
  out$downstream <- downstream
  out$bin_size <- bin_size
  
  out
}
### get the folder sample ###
sample_folder = list.files("../All_protein_coding_TSS_output",full.names = T)


# combine the files 
for (folder in sample_folder) {
  
  #print(folder)
  ### get the files
  files <- list.files(
    path = folder,
    pattern = "\\.tab$",
    full.names = TRUE)
  
  # get the chromosome location
  chr = sub("^chr", "", sub("_.*", "", basename(folder)))
  
  # get the loss sample
  
  # get loss chromosome 
  chr_loss_df = TACG_df[TACG_df[,chr] == -1 &TACG_df$Genome_doublings == 0 ,]
  
  chr_loss_sample = chr_loss_df$Sample
  chr_loss_sample = chr_loss_sample[!is.na(chr_loss_sample)]
  
  # to get the ATACseq_sample
  sample_atac %>% 
    dplyr::mutate(sample_id = str_sub(aliquot_id, 1,15)) %>% 
    dplyr::filter(sample_id %in% chr_loss_sample) %>% 
    dplyr::mutate(bam_prefix = gsub("-", "_", bam_prefix)) %>% 
    dplyr::select(bam_prefix) %>% 
    unlist(use.names = F) -> chr_loss_atac_sample
  
  # get the loss files
  files_chr_loss_selected <- files[basename(files) %in% paste(chr_loss_atac_sample, ".tab", sep = "")]
  
  # read the loss files 
  chr_loss_profiles <- bind_rows(lapply(files_chr_loss_selected, read_profile_tss))
  
  # ----------------------------
  # average across loss files
  # ----------------------------
  summary_loss_df <- chr_loss_profiles %>%
    group_by(group, bin) %>%
    summarise(
      mean = mean(signal, na.rm = TRUE),
      sd   = sd(signal, na.rm = TRUE),
      n    = n(),
      se   = sd / sqrt(n),
      .groups = "drop"
    )
  
  
  # get the neutral files 
  # get loss chromosome 
  chr_neutral_df = TACG_df[TACG_df[,chr] == 0 &TACG_df$Genome_doublings == 0 ,]
  
  chr_neutral_sample = chr_neutral_df$Sample
  chr_neutral_sample = chr_neutral_sample[!is.na(chr_neutral_sample)]
  
  # to get the ATACseq_sample
  sample_atac %>% 
    dplyr::mutate(sample_id = str_sub(aliquot_id, 1,15)) %>% 
    dplyr::filter(sample_id %in% chr_neutral_sample) %>% 
    dplyr::mutate(bam_prefix = gsub("-", "_", bam_prefix)) %>% 
    dplyr::select(bam_prefix) %>% 
    unlist(use.names = F) -> chr_neutral_atac_sample
  
  # get the loss files
  files_chr_neutral_selected <- files[basename(files) %in% paste(chr_neutral_atac_sample, ".tab", sep = "")]
  
  # read the loss files 
  chr_neutral_profiles <- bind_rows(lapply(files_chr_neutral_selected, read_profile_tss))
  
  
  # get all the neutral files 
  summary_neutral_df <- chr_neutral_profiles %>%
    group_by(group, bin) %>%
    summarise(
      mean = mean(signal, na.rm = TRUE),
      sd   = sd(signal, na.rm = TRUE),
      n    = n(),
      se   = sd / sqrt(n),
      .groups = "drop"
    )
  
  
  # normalise the copy number
  summary_neutral_df %>% 
    dplyr::mutate(mean = mean/ 2) -> summary_neutral_df_mdf
  
  summary_neutral_df_mdf$condition = "neutral"
  summary_loss_df$condition = "loss"
  
  # combine the table 
  df_cmb = rbind(summary_neutral_df_mdf,summary_loss_df)
  
  head(df_cmb)
  df_cmb$chr_position = paste("chr", chr, sep = "")
  # ----------------------------
  # axis info from first file
  # ----------------------------
  
  # save the file into big formation
  if(!exists("comb_df")){
    comb_df = df_cmb
  }else{
    comb_df = rbind(comb_df, df_cmb)
  }
  
  
}


h2 <- sub("^#", "", readLines(files_chr_loss_selected[1], n = 2)[2])

upstream   <- as.numeric(sub(".*upstream:(\\d+).*", "\\1", h2))
downstream <- as.numeric(sub(".*downstream:(\\d+).*", "\\1", h2))
bin_size   <- as.numeric(sub(".*bin size:(\\d+).*", "\\1", h2))

up_bins   <- upstream / bin_size
down_bins <- downstream / bin_size
total_bins <- up_bins + down_bins

tss_pos <- up_bins + 1


plot_df <- comb_df %>%
  #dplyr::filter(chr == "chr8p")
  dplyr::group_by(bin, condition) %>%
  dplyr::summarise(
    median_value = median(mean, na.rm = TRUE),
    n_files      = sum(!is.na(mean)),
    se_value     = sd(mean, na.rm = TRUE) / sqrt(n_files),
    .groups = "drop"
  )


# plot the figure
library(ggthemes)
ggplot(
  plot_df,
  aes(
    x = bin,
    y = median_value,
    color = condition,
    fill = condition
  )
) +
  geom_ribbon(
    aes(
      ymin = median_value - se_value,
      ymax = median_value + se_value
    ),
    color = NA,
    alpha = 0.20
  ) +
  geom_smooth(
    method = "loess",
    span = 0.05,
    se = FALSE,
    linewidth = 1.1
  )+
  geom_vline(
    xintercept = tss_pos,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(values = c(
    "loss"    = "#FCAE1E",
    "neutral" = "#7852A9"
  )) +
  scale_fill_manual(values = c(
    "loss"    = "#FCAE1E",
    "neutral" = "#7852A9"
  )) +
  scale_x_continuous(
    breaks = c(1, tss_pos, total_bins),
    labels = c(
      paste0("-", upstream / 1000, " Kb"),
      "TSS",
      paste0("+", downstream / 1000, " Kb")
    )
  ) +
  #coord_cartesian(ylim = c(0, 25)) +
  labs(
    x = NULL,
    y = "Signal",
    color = NULL,
    fill = NULL
  ) +
  geom_rangeframe(color = "black") +
  theme_tufte(base_size = 12)+
  #theme_classic(base_size = 18) +
  theme(
    #  legend.position = c(0.5, 0.85),
    #axis.line = element_line(color =  "black"),
    #  axis.ticks = element_line(linewidth = 0.8),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) 

############################################
#################Figure 3C #################
############################################
# library
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)

##data
sample_atac = read.delim("../TCGA_identifier_mapping.txt")
TACG_df = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Aneuploidy patient datasets/mmc2.xlsx", skip = 1)


# function 
# ----------------------------
# function to read one matrix
# ----------------------------
read_profile_tss <- function(file) {
  
  x <- readLines(file)
  
  # line 1: group names and group sizes
  h1 <- sub("^#", "", x[1])
  ginfo <- strsplit(h1, "\t")[[1]]
  
  gnames <- sub(":(.*)$", "", ginfo)
  gsizes <- as.numeric(sub("^.*:", "", ginfo))
  
  # line 2: parameters
  h2 <- sub("^#", "", x[2])
  
  upstream   <- as.numeric(sub(".*upstream:(\\d+).*", "\\1", h2))
  downstream <- as.numeric(sub(".*downstream:(\\d+).*", "\\1", h2))
  bin_size   <- as.numeric(sub(".*bin size:(\\d+).*", "\\1", h2))
  
  # numeric matrix starts after first 3 lines
  num_lines <- x[-c(1, 2, 3)]
  
  mat <- do.call(
    rbind,
    lapply(num_lines, function(z) as.numeric(strsplit(z, "\t")[[1]]))
  )
  
  mat <- as.matrix(mat)
  
  group <- rep(gnames, gsizes)
  
  out <- do.call(
    rbind,
    lapply(gnames, function(g) {
      idx <- which(group == g)
      
      data.frame(
        group = g,
        bin = seq_len(ncol(mat)),
        signal = colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
      )
    })
  )
  
  out$sample <- basename(file)
  out$upstream <- upstream
  out$downstream <- downstream
  out$bin_size <- bin_size
  
  out
}
### get the folder sample ###
sample_folder = list.files("../Outputs_data",full.names = T)

for (folder in sample_folder) {
  
  ### get the files
  files <- list.files(
    path = folder,
    pattern = "\\.tab$",
    full.names = TRUE)
  
  # get the chromosome location
  chr = str_sub(folder, start = 27, end = nchar(folder))
  
  # get the loss sample
  
  # get loss chromosome 
  chr_loss_df = TACG_df[TACG_df[,chr] == -1 &TACG_df$Genome_doublings == 0 ,]
  
  chr_loss_sample = chr_loss_df$Sample
  chr_loss_sample = chr_loss_sample[!is.na(chr_loss_sample)]
  
  # to get the ATACseq_sample
  sample_atac %>% 
    dplyr::mutate(sample_id = str_sub(aliquot_id, 1,15)) %>% 
    dplyr::filter(sample_id %in% chr_loss_sample) %>% 
    dplyr::mutate(bam_prefix = gsub("-", "_", bam_prefix)) %>% 
    dplyr::select(bam_prefix) %>% 
    unlist(use.names = F) -> chr_loss_atac_sample
  
  # get the loss files
  files_chr_loss_selected <- files[basename(files) %in% paste(chr_loss_atac_sample, ".tab", sep = "")]
  
  # read the loss files 
  chr_loss_profiles <- bind_rows(lapply(files_chr_loss_selected, read_profile_tss))
  
  # ----------------------------
  # average across loss files
  # ----------------------------
  summary_loss_df <- chr_loss_profiles %>%
    group_by(group, bin) %>%
    summarise(
      mean = mean(signal, na.rm = TRUE),
      sd   = sd(signal, na.rm = TRUE),
      n    = n(),
      se   = sd / sqrt(n),
      .groups = "drop"
    )
  
  
  # get the neutral files 
  # get loss chromosome 
  chr_neutral_df = TACG_df[TACG_df[,chr] == 0 &TACG_df$Genome_doublings == 0 ,]
  
  chr_neutral_sample = chr_neutral_df$Sample
  chr_neutral_sample = chr_neutral_sample[!is.na(chr_neutral_sample)]
  
  # to get the ATACseq_sample
  sample_atac %>% 
    dplyr::mutate(sample_id = str_sub(aliquot_id, 1,15)) %>% 
    dplyr::filter(sample_id %in% chr_neutral_sample) %>% 
    dplyr::mutate(bam_prefix = gsub("-", "_", bam_prefix)) %>% 
    dplyr::select(bam_prefix) %>% 
    unlist(use.names = F) -> chr_neutral_atac_sample
  
  # get the loss files
  files_chr_neutral_selected <- files[basename(files) %in% paste(chr_neutral_atac_sample, ".tab", sep = "")]
  
  # read the loss files 
  chr_neutral_profiles <- bind_rows(lapply(files_chr_neutral_selected, read_profile_tss))
  
  
  # get all the neutral files 
  summary_neutral_df <- chr_neutral_profiles %>%
    group_by(group, bin) %>%
    summarise(
      mean = mean(signal, na.rm = TRUE),
      sd   = sd(signal, na.rm = TRUE),
      n    = n(),
      se   = sd / sqrt(n),
      .groups = "drop"
    )
  
  
  # normalise the copy number
  summary_neutral_df %>% 
    dplyr::mutate(mean = mean/ 2) -> summary_neutral_df_mdf
  
  summary_neutral_df_mdf$condition = "neutral"
  summary_loss_df$condition = "loss"
  
  # combine the table 
  df_cmb = rbind(summary_neutral_df_mdf,summary_loss_df)
  
  # 8p file is different 
  if(chr == "8p"){
    df_cmb$chr = "chr8p"
    df_cmb %>% 
      mutate(gene_label = ifelse(group == "PPP2CB.bed", "sig", "nonsig")) ->df_cmb
  }else{
    # add the gene label
    df_cmb$gene_label = sub("^.*_(.*?)\\.bed$", "\\1", df_cmb$group)
    
    # add chromosome information
    df_cmb$chr = sub("_.*$", "", df_cmb$group)
    # check the tbl
    head(df_cmb)
  }
  
  
  # save the file into big formation
  if(!exists("comb_df")){
    comb_df = df_cmb
  }else{
    comb_df = rbind(comb_df, df_cmb)
  }
  
}


h2 <- sub("^#", "", readLines(files_chr_loss_selected[1], n = 2)[2])

upstream   <- as.numeric(sub(".*upstream:(\\d+).*", "\\1", h2))
downstream <- as.numeric(sub(".*downstream:(\\d+).*", "\\1", h2))
bin_size   <- as.numeric(sub(".*bin size:(\\d+).*", "\\1", h2))

up_bins   <- upstream / bin_size
down_bins <- downstream / bin_size
total_bins <- up_bins + down_bins

tss_pos <- up_bins + 1

#comb_df_sig <- comb_df %>%
# dplyr::filter(gene_label == "sig")

plot_df <- comb_df %>%
  #dplyr::filter(chr == "chr8p")
  dplyr::group_by(gene_label, bin, condition) %>%
  dplyr::summarise(
    median_value = median(mean, na.rm = TRUE),
    n_files      = sum(!is.na(mean)),
    se_value     = sd(mean, na.rm = TRUE) / sqrt(n_files),
    .groups = "drop"
  )


# plot the figure
library(ggthemes)
ggplot(
  plot_df,
  aes(
    x = bin,
    y = median_value,
    color = condition,
    fill = condition
  )
) +
  geom_ribbon(
    aes(
      ymin = median_value - se_value,
      ymax = median_value + se_value
    ),
    color = NA,
    alpha = 0.20
  ) +
  geom_smooth(
    method = "loess",
    span = 0.05,
    se = FALSE,
    linewidth = 1.1
  )+
  geom_vline(
    xintercept = tss_pos,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(values = c(
    "loss"    = "#FCAE1E",
    "neutral" = "#7852A9"
  )) +
  scale_fill_manual(values = c(
    "loss"    = "#FCAE1E",
    "neutral" = "#7852A9"
  )) +
  scale_x_continuous(
    breaks = c(1, tss_pos, total_bins),
    labels = c(
      paste0("-", upstream / 1000, " Kb"),
      "TSS",
      paste0("+", downstream / 1000, " Kb")
    )
  ) +
  #coord_cartesian(ylim = c(0, 25)) +
  labs(
    x = NULL,
    y = "Signal",
    color = NULL,
    fill = NULL
  ) +
  geom_rangeframe(color = "black") +
  theme_tufte(base_size = 12)+
  #theme_classic(base_size = 18) +
  theme(
    #  legend.position = c(0.5, 0.85),
    #axis.line = element_line(color =  "black"),
    #  axis.ticks = element_line(linewidth = 0.8),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) +
  facet_wrap(~gene_label)


# 800, 300
############################################
#################Figure 3D #################
############################################

####### look at the inidvidual chromosome arm loss ##########
# Chr8p loss#
plot_df_8p <- comb_df %>%
  filter(chr == "chr8p", is.finite(bin), is.finite(mean)) %>%
  group_by(gene_label, condition) %>%
  group_modify(~ {
    df <- arrange(.x, bin)
    
    fit <- loess(
      mean ~ bin,
      data = df,
      span = 0.05,
      na.action = na.exclude
    )
    
    df %>%
      mutate(
        mean_smooth = as.numeric(
          predict(fit, newdata = data.frame(bin = bin))
        ),
        lower = mean_smooth - 1.5 * se,
        upper = mean_smooth + 1.5 * se
      )
  }) %>%
  ungroup()

ggplot(
  plot_df_8p,
  aes(
    x = bin,
    y = mean_smooth,
    color = condition,
    fill = condition,
    group = condition
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    color = NA,
    alpha = 0.20
  ) +
  geom_line(linewidth = 1.1) +
  geom_vline(
    xintercept = tss_pos,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(values = c(
    loss = "#FCAE1E",
    neutral = "#7852A9"
  )) +
  scale_fill_manual(values = c(
    loss = "#FCAE1E",
    neutral = "#7852A9"
  )) +
  scale_x_continuous(
    breaks = c(1, tss_pos, total_bins),
    labels = c(
      paste0("-", upstream / 1000, " Kb"),
      "TSS",
      paste0("+", downstream / 1000, " Kb")
    )
  ) +
  labs(
    x = NULL,
    y = "Signal",
    color = NULL,
    fill = NULL
  ) +
  ggthemes::geom_rangeframe(color = "black") +
  ggthemes::theme_tufte(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  facet_wrap(~gene_label)


########### 5q ###########
plot_df_5q <- comb_df %>%
  filter(chr == "chr5q", is.finite(bin), is.finite(mean)) %>%
  group_by(gene_label, condition) %>%
  group_modify(~ {
    df <- arrange(.x, bin)
    
    fit <- loess(
      mean ~ bin,
      data = df,
      span = 0.05,
      na.action = na.exclude
    )
    
    df %>%
      mutate(
        mean_smooth = as.numeric(
          predict(fit, newdata = data.frame(bin = bin))
        ),
        lower = mean_smooth - 1.5 * se,
        upper = mean_smooth + 1.5 * se
      )
  }) %>%
  ungroup()

ggplot(
  plot_df_5q,
  aes(
    x = bin,
    y = mean_smooth,
    color = condition,
    fill = condition,
    group = condition
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    color = NA,
    alpha = 0.20
  ) +
  geom_line(linewidth = 1.1) +
  geom_vline(
    xintercept = tss_pos,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(values = c(
    loss = "#FCAE1E",
    neutral = "#7852A9"
  )) +
  scale_fill_manual(values = c(
    loss = "#FCAE1E",
    neutral = "#7852A9"
  )) +
  scale_x_continuous(
    breaks = c(1, tss_pos, total_bins),
    labels = c(
      paste0("-", upstream / 1000, " Kb"),
      "TSS",
      paste0("+", downstream / 1000, " Kb")
    )
  ) +
  labs(
    x = NULL,
    y = "Signal",
    color = NULL,
    fill = NULL
  ) +
  ggthemes::geom_rangeframe(color = "black") +
  ggthemes::theme_tufte(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  facet_wrap(~gene_label)




############################################
#################Figure 3E #################
############################################
# Project: H3K27ac HiChIP Peaks Across Recurrent Chromosome-Loss Events in Cancer for all protein coding genes

# load the library
library(tidyverse)
library(ggthemes)
library(ggpubr)

# laod the data from the HiChIP
# The original paper is from "Three-dimensional genome landscape of primary human cancers, 2025"

# Load the data
# Metasheet data for the HiChIP
# get all the sample info 

sample_meta_HiChIP = read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/HiChIP dataset/TCGA_HiChIP_metadata.txt")

# TCGA sample aneuploidy
TACG_df = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Aneuploidy patient datasets/mmc2.xlsx", skip = 1)


# function 
# ----------------------------
# function to read one matrix
# ----------------------------
read_profile_tss <- function(file) {
  
  x <- readLines(file)
  
  # line 1: group names and group sizes
  h1 <- sub("^#", "", x[1])
  ginfo <- strsplit(h1, "\t")[[1]]
  
  gnames <- sub(":(.*)$", "", ginfo)
  gsizes <- as.numeric(sub("^.*:", "", ginfo))
  
  # line 2: parameters
  h2 <- sub("^#", "", x[2])
  
  upstream   <- as.numeric(sub(".*upstream:(\\d+).*", "\\1", h2))
  downstream <- as.numeric(sub(".*downstream:(\\d+).*", "\\1", h2))
  bin_size   <- as.numeric(sub(".*bin size:(\\d+).*", "\\1", h2))
  
  # numeric matrix starts after first 3 lines
  num_lines <- x[-c(1, 2, 3)]
  
  mat <- do.call(
    rbind,
    lapply(num_lines, function(z) as.numeric(strsplit(z, "\t")[[1]]))
  )
  
  mat <- as.matrix(mat)
  
  group <- rep(gnames, gsizes)
  
  out <- do.call(
    rbind,
    lapply(gnames, function(g) {
      idx <- which(group == g)
      
      data.frame(
        group = g,
        bin = seq_len(ncol(mat)),
        signal = colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
      )
    })
  )
  
  out$sample <- basename(file)
  out$upstream <- upstream
  out$downstream <- downstream
  out$bin_size <- bin_size
  
  out
}
### get the folder sample ###
sample_folder = list.files("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/HiChIP dataset/All_protein_coding_TSS_output",full.names = T)



for (folder in sample_folder) {
  
  ### get the files
  files <- list.files(
    path = folder,
    pattern = "\\.tab$",
    full.names = TRUE)
  
  # get the chromosome location
  chr <- str_extract(folder, "(?:[1-9]|1[0-9]|2[0-2]|X|Y)[pq]")
  
  if(is.na(chr)){
    chr = "17p"
  }
  
  
  print(chr)
  # get the loss sample
  
  # get loss chromosome 
  chr_loss_df = TACG_df[TACG_df[,chr] == -1 &TACG_df$Genome_doublings == 0 ,]
  
  chr_loss_sample = chr_loss_df$Sample
  chr_loss_sample = chr_loss_sample[!is.na(chr_loss_sample)]
  
  # to get the ATACseq_sample
  sample_meta_HiChIP %>% 
    dplyr::mutate(sample_id = str_sub(aliquot_submitter_id, 1,15)) %>% 
    dplyr::filter(sample_id %in% chr_loss_sample) %>% 
    dplyr::mutate(bam_prefix = sub(
      "(_H3K27ac).*$",
      "\\1",
      basename(file_name)
    )) %>% 
    dplyr::select(bam_prefix) %>% 
    distinct(bam_prefix) %>% 
    unlist(use.names = F) -> chr_loss_hiChip_sample
  
  if(length(chr_loss_hiChip_sample) >= 5){
    # get the loss files
    files_chr_loss_selected <- files[basename(files) %in% paste(chr_loss_hiChip_sample, "_1D-signal.tab", sep = "")]
    
    # read the loss files 
    chr_loss_profiles <- bind_rows(lapply(files_chr_loss_selected, read_profile_tss))
    
    # ----------------------------
    # average across loss files
    # ----------------------------
    summary_loss_df <- chr_loss_profiles %>%
      group_by(group, bin) %>%
      summarise(
        mean = mean(signal, na.rm = TRUE),
        sd   = sd(signal, na.rm = TRUE),
        n    = n(),
        se   = sd / sqrt(n),
        .groups = "drop"
      )
    
    
    # get the neutral files 
    # get loss chromosome 
    chr_neutral_df = TACG_df[TACG_df[,chr] == 0 &TACG_df$Genome_doublings == 0 ,]
    
    chr_neutral_sample = chr_neutral_df$Sample
    chr_neutral_sample = chr_neutral_sample[!is.na(chr_neutral_sample)]
    
    # to get the ATACseq_sample
    sample_meta_HiChIP %>% 
      dplyr::mutate(sample_id = str_sub(aliquot_submitter_id, 1,15)) %>% 
      dplyr::filter(sample_id %in% chr_neutral_sample) %>% 
      dplyr::mutate(bam_prefix = sub(
        "(_H3K27ac).*$",
        "\\1",
        basename(file_name)
      )) %>% 
      dplyr::select(bam_prefix) %>% 
      distinct(bam_prefix) %>% 
      unlist(use.names = F) -> chr_neutral_hiChip_sample
    
    # get the loss files
    files_chr_neutral_selected <- files[basename(files) %in% paste(chr_neutral_hiChip_sample, "_1D-signal.tab", sep = "")]
    
    # read the loss files 
    chr_neutral_profiles <- bind_rows(lapply(files_chr_neutral_selected, read_profile_tss))
    
    
    # get all the neutral files 
    summary_neutral_df <- chr_neutral_profiles %>%
      group_by(group, bin) %>%
      summarise(
        mean = mean(signal, na.rm = TRUE),
        sd   = sd(signal, na.rm = TRUE),
        n    = n(),
        se   = sd / sqrt(n),
        .groups = "drop"
      )
    
    
    # normalise the copy number
    summary_neutral_df %>% 
      dplyr::mutate(mean = mean/ 2) -> summary_neutral_df_mdf
    
    summary_neutral_df_mdf$condition = "neutral"
    summary_loss_df$condition = "loss"
    
    # combine the table 
    df_cmb = rbind(summary_neutral_df_mdf,summary_loss_df)
    
    
    
    # add the gene label
    df_cmb$gene_label = sub("^.*_(.*?)\\.bed$", "\\1", df_cmb$group)
    
    # add chromosome information
    df_cmb$chr = chr
    # check the tbl
    head(df_cmb)
    
    # save the file into big formation
    if(!exists("comb_df")){
      comb_df = df_cmb
    }else{
      comb_df = rbind(comb_df, df_cmb)
    }
    
  }
  
  
}


h2 <- sub("^#", "", readLines(files_chr_loss_selected[1], n = 2)[2])

upstream   <- as.numeric(sub(".*upstream:(\\d+).*", "\\1", h2))
downstream <- as.numeric(sub(".*downstream:(\\d+).*", "\\1", h2))
bin_size   <- as.numeric(sub(".*bin size:(\\d+).*", "\\1", h2))

up_bins   <- upstream / bin_size
down_bins <- downstream / bin_size
total_bins <- up_bins + down_bins

tss_pos <- up_bins + 1

#comb_df_sig <- comb_df %>%
# dplyr::filter(gene_label == "chr17")


#The missing ribbon occurs because sd() returns NA when a bin has only one observation.
#Also, completely absent bins do not appear in plot_df. 
#Complete all bins first, then interpolate missing median and SE values within each gene_label × condition group.
# Interpolate missing values
interpolate_na <- function(x, y) {
  valid <- is.finite(x) & is.finite(y)
  
  if (sum(valid) >= 2) {
    approx(
      x = x[valid],
      y = y[valid],
      xout = x,
      rule = 2
    )$y
  } else if (sum(valid) == 1) {
    rep(y[valid][1], length(y))
  } else {
    rep(NA_real_, length(y))
  }
}

# Calculate mean and SE for each bin, then fill missing values
plot_df <- comb_df %>%
  dplyr::group_by(condition, bin) %>%
  dplyr::summarise(
    mean_value = mean(mean, na.rm = TRUE),
    n_files = sum(is.finite(mean)),
    se = {
      values <- mean[is.finite(mean)]
      
      if (length(values) >= 2) {
        sd(values) / sqrt(length(values))
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    mean_value = dplyr::na_if(mean_value, NaN)
  ) %>%
  dplyr::group_by( condition) %>%
  tidyr::complete(bin = seq_len(total_bins)) %>%
  dplyr::arrange(bin, .by_group = TRUE) %>%
  dplyr::mutate(
    mean_value = interpolate_na(bin, mean_value),
    se = interpolate_na(bin, se),
    ribbon_lower = pmax(0, mean_value - se),
    ribbon_upper = mean_value + se
  ) %>%
  dplyr::ungroup()

ggplot(
  plot_df,
  aes(
    x = bin,
    y = mean_value,
    color = condition,
    fill = condition
  )
) +
  geom_ribbon(
    aes(
      ymin = ribbon_lower,
      ymax = ribbon_upper,
      group = condition
    ),
    color = NA,
    alpha = 0.20
  ) +
  geom_smooth(
    aes(group = condition),
    method = "loess",
    span = 0.05,
    se = FALSE,
    linewidth = 1.1
  ) +
  geom_vline(
    xintercept = tss_pos,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = c(
      "loss" = "#FCAE1E",
      "neutral" = "#7852A9"
    )
  ) +
  scale_fill_manual(
    values = c(
      "loss" = "#FCAE1E",
      "neutral" = "#7852A9"
    )
  ) +
  scale_x_continuous(
    breaks = c(1, tss_pos, total_bins),
    labels = c(
      paste0("-", upstream / 1000, " Kb"),
      "TSS",
      paste0("+", downstream / 1000, " Kb")
    )
  ) +
  labs(
    x = NULL,
    y = "Signal",
    color = NULL,
    fill = NULL
  ) +
  geom_rangeframe(color = "black") +
  theme_tufte(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) +
  scale_y_continuous(limits = c(0, 3),
                     breaks = seq(0, 3, by = 1))



############################################
#################Figure 3F #################
############################################
# Project: H3K27ac HiChIP Peaks Across Recurrent Chromosome-Loss Events in Cancer

# load the library
library(tidyverse)
library(ggthemes)
library(ggpubr)

# laod the data from the HiChIP
# The original paper is from "Three-dimensional genome landscape of primary human cancers, 2025"

# Load the data
# Metasheet data for the HiChIP
# get all the sample info 

sample_meta_HiChIP = read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/HiChIP dataset/TCGA_HiChIP_metadata.txt")

# TCGA sample aneuploidy
TACG_df = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Aneuploidy patient datasets/mmc2.xlsx", skip = 1)


# function 
# ----------------------------
# function to read one matrix
# ----------------------------
read_profile_tss <- function(file) {
  
  x <- readLines(file)
  
  # line 1: group names and group sizes
  h1 <- sub("^#", "", x[1])
  ginfo <- strsplit(h1, "\t")[[1]]
  
  gnames <- sub(":(.*)$", "", ginfo)
  gsizes <- as.numeric(sub("^.*:", "", ginfo))
  
  # line 2: parameters
  h2 <- sub("^#", "", x[2])
  
  upstream   <- as.numeric(sub(".*upstream:(\\d+).*", "\\1", h2))
  downstream <- as.numeric(sub(".*downstream:(\\d+).*", "\\1", h2))
  bin_size   <- as.numeric(sub(".*bin size:(\\d+).*", "\\1", h2))
  
  # numeric matrix starts after first 3 lines
  num_lines <- x[-c(1, 2, 3)]
  
  mat <- do.call(
    rbind,
    lapply(num_lines, function(z) as.numeric(strsplit(z, "\t")[[1]]))
  )
  
  mat <- as.matrix(mat)
  
  group <- rep(gnames, gsizes)
  
  out <- do.call(
    rbind,
    lapply(gnames, function(g) {
      idx <- which(group == g)
      
      data.frame(
        group = g,
        bin = seq_len(ncol(mat)),
        signal = colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
      )
    })
  )
  
  out$sample <- basename(file)
  out$upstream <- upstream
  out$downstream <- downstream
  out$bin_size <- bin_size
  
  out
}
### get the folder sample ###
sample_folder = list.files("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/HiChIP dataset/H3K27ac_all_chromosome",full.names = T)



for (folder in sample_folder) {
  
  ### get the files
  files <- list.files(
    path = folder,
    pattern = "\\.tab$",
    full.names = TRUE)
  
  # get the chromosome location
  chr <- str_extract(folder, "(?:[1-9]|1[0-9]|2[0-2]|X|Y)[pq]")
  
  if(is.na(chr)){
    chr = "17p"
  }
  
  
  print(chr)
  # get the loss sample
  
  # get loss chromosome 
  chr_loss_df = TACG_df[TACG_df[,chr] == -1 &TACG_df$Genome_doublings == 0 ,]
  
  chr_loss_sample = chr_loss_df$Sample
  chr_loss_sample = chr_loss_sample[!is.na(chr_loss_sample)]
  
  # to get the ATACseq_sample
  sample_meta_HiChIP %>% 
    dplyr::mutate(sample_id = str_sub(aliquot_submitter_id, 1,15)) %>% 
    dplyr::filter(sample_id %in% chr_loss_sample) %>% 
    dplyr::mutate(bam_prefix = sub(
      "(_H3K27ac).*$",
      "\\1",
      basename(file_name)
    )) %>% 
    dplyr::select(bam_prefix) %>% 
    distinct(bam_prefix) %>% 
    unlist(use.names = F) -> chr_loss_hiChip_sample
  
  if(length(chr_loss_hiChip_sample) >= 5){
    # get the loss files
    files_chr_loss_selected <- files[basename(files) %in% paste(chr_loss_hiChip_sample, "_1D-signal.tab", sep = "")]
    
    # read the loss files 
    chr_loss_profiles <- bind_rows(lapply(files_chr_loss_selected, read_profile_tss))
    
    # ----------------------------
    # average across loss files
    # ----------------------------
    summary_loss_df <- chr_loss_profiles %>%
      group_by(group, bin) %>%
      summarise(
        mean = mean(signal, na.rm = TRUE),
        sd   = sd(signal, na.rm = TRUE),
        n    = n(),
        se   = sd / sqrt(n),
        .groups = "drop"
      )
    
    
    # get the neutral files 
    # get loss chromosome 
    chr_neutral_df = TACG_df[TACG_df[,chr] == 0 &TACG_df$Genome_doublings == 0 ,]
    
    chr_neutral_sample = chr_neutral_df$Sample
    chr_neutral_sample = chr_neutral_sample[!is.na(chr_neutral_sample)]
    
    # to get the ATACseq_sample
    sample_meta_HiChIP %>% 
      dplyr::mutate(sample_id = str_sub(aliquot_submitter_id, 1,15)) %>% 
      dplyr::filter(sample_id %in% chr_neutral_sample) %>% 
      dplyr::mutate(bam_prefix = sub(
        "(_H3K27ac).*$",
        "\\1",
        basename(file_name)
      )) %>% 
      dplyr::select(bam_prefix) %>% 
      distinct(bam_prefix) %>% 
      unlist(use.names = F) -> chr_neutral_hiChip_sample
    
    # get the loss files
    files_chr_neutral_selected <- files[basename(files) %in% paste(chr_neutral_hiChip_sample, "_1D-signal.tab", sep = "")]
    
    # read the loss files 
    chr_neutral_profiles <- bind_rows(lapply(files_chr_neutral_selected, read_profile_tss))
    
    
    # get all the neutral files 
    summary_neutral_df <- chr_neutral_profiles %>%
      group_by(group, bin) %>%
      summarise(
        mean = mean(signal, na.rm = TRUE),
        sd   = sd(signal, na.rm = TRUE),
        n    = n(),
        se   = sd / sqrt(n),
        .groups = "drop"
      )
    
    
    # normalise the copy number
    summary_neutral_df %>% 
      dplyr::mutate(mean = mean/ 2) -> summary_neutral_df_mdf
    
    summary_neutral_df_mdf$condition = "neutral"
    summary_loss_df$condition = "loss"
    
    # combine the table 
    df_cmb = rbind(summary_neutral_df_mdf,summary_loss_df)
    
    
    
    # add the gene label
    df_cmb$gene_label = sub("^.*_(.*?)\\.bed$", "\\1", df_cmb$group)
    
    # add chromosome information
    df_cmb$chr = sub("_.*$", "", df_cmb$group)
    # check the tbl
    head(df_cmb)
    
    # save the file into big formation
    if(!exists("comb_df")){
      comb_df = df_cmb
    }else{
      comb_df = rbind(comb_df, df_cmb)
    }
    
  }
  
  
}



h2 <- sub("^#", "", readLines(files_chr_loss_selected[1], n = 2)[2])

upstream   <- as.numeric(sub(".*upstream:(\\d+).*", "\\1", h2))
downstream <- as.numeric(sub(".*downstream:(\\d+).*", "\\1", h2))
bin_size   <- as.numeric(sub(".*bin size:(\\d+).*", "\\1", h2))

up_bins   <- upstream / bin_size
down_bins <- downstream / bin_size
total_bins <- up_bins + down_bins

tss_pos <- up_bins + 1

#comb_df_sig <- comb_df %>%
# dplyr::filter(gene_label == "chr17")


#The missing ribbon occurs because sd() returns NA when a bin has only one observation.
#Also, completely absent bins do not appear in plot_df. 
#Complete all bins first, then interpolate missing median and SE values within each gene_label × condition group.
interpolate_na <- function(x, y) {
  valid <- is.finite(y)
  
  if (sum(valid) >= 2) {
    approx(
      x = x[valid],
      y = y[valid],
      xout = x,
      rule = 2
    )$y
  } else if (sum(valid) == 1) {
    rep(y[valid][1], length(y))
  } else {
    rep(NA_real_, length(y))
  }
}

plot_df <- comb_df %>%
  dplyr::group_by(gene_label, bin, condition) %>%
  dplyr::summarise(
    median_value = median(mean, na.rm = TRUE),
    n_files = sum(!is.na(mean)),
    
    # SE cannot be calculated from only one sample
    se_value = dplyr::if_else(
      n_files >= 2,
      sd(mean, na.rm = TRUE) / sqrt(n_files),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    median_value = dplyr::na_if(median_value, NaN)
  ) %>%
  
  # Add bins that are completely absent from the original data
  dplyr::group_by(gene_label, condition) %>%
  tidyr::complete(bin = seq_len(total_bins)) %>%
  dplyr::arrange(bin, .by_group = TRUE) %>%
  
  # Interpolate missing values within each curve
  dplyr::mutate(
    median_value = interpolate_na(bin, median_value),
    se_value = interpolate_na(bin, se_value),
    ribbon_lower = pmax(0, median_value - se_value),
    ribbon_upper = median_value + se_value
  ) %>%
  dplyr::ungroup()



ggplot(
  plot_df,
  aes(
    x = bin,
    y = median_value,
    color = condition,
    fill = condition
  )
) +
  geom_ribbon(
    aes(
      ymin = ribbon_lower,
      ymax = ribbon_upper,
      group = condition
    ),
    color = NA,
    alpha = 0.20
  ) +
  geom_smooth(
    aes(group = condition),
    method = "loess",
    span = 0.05,
    se = FALSE,
    linewidth = 1.1
  ) +
  geom_vline(
    xintercept = tss_pos,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = c(
      "loss" = "#FCAE1E",
      "neutral" = "#7852A9"
    )
  ) +
  scale_fill_manual(
    values = c(
      "loss" = "#FCAE1E",
      "neutral" = "#7852A9"
    )
  ) +
  scale_x_continuous(
    breaks = c(1, tss_pos, total_bins),
    labels = c(
      paste0("-", upstream / 1000, " Kb"),
      "TSS",
      paste0("+", downstream / 1000, " Kb")
    )
  ) +
  labs(
    x = NULL,
    y = "Signal",
    color = NULL,
    fill = NULL
  ) +
  geom_rangeframe(color = "black") +
  theme_tufte(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) +
  facet_wrap(~gene_label)

# 800, 300

############################################
#################Figure 3G #################
############################################
############################################################################################################################################################
####  in 37 sig: #Do significant paralogs on lost chromosomes exhibit less dosage compensation compared to non-significant paralogs on lost chromosomes?####
############################################################################################################################################################
sig_37_paralog_df %>% 
  
  mutate(
    # pull captures per row
    m = str_match(aneuploid_loss_chr, "^chr(\\d+|X|Y)_([pq])$"),
    repl = paste0("X", m[,2], m[,3]),
    # keep chrY as-is; otherwise use replacement
    aneuploid_loss_chr = if_else(aneuploid_loss_chr == "chrY", aneuploid_loss_chr, repl)
  ) %>% 
  dplyr::select(paralog_gene,aneuploid_loss_chr) -> sig_loc_gene_df

# non sig loc gene and pos
non_sig_para_query_tbl %>% 
  dplyr::select(para_gene_1,chr_position) -> non_sig_query_1
colnames(non_sig_query_1) = c("gene_names", "chr_pos")
non_sig_para_query_tbl %>% 
  dplyr::select(para_gene_2,chr_para2_position) -> non_sig_query_2
colnames(non_sig_query_2) = c("gene_names", "chr_pos")

non_sig_query_list = rbind(non_sig_query_1, non_sig_query_2)
# can not include chrX loss
non_sig_query_list %>% 
  distinct(gene_names, .keep_all = T) %>% 
  dplyr::filter(!(gene_names %in% sig_loc_gene_df$paralog_gene)) %>% 
  dplyr::filter(chr_pos != "chrX") %>% 
  dplyr::filter(gene_names %in% colnames(CCLE_exp))-> non_sig_query_list_flt

# CCLE cell line with each chromsome arm gain or loss info
#df_aneu_score = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/arm_call_scores.csv")
#get all the aneuploidy cells
df_aneu_score %>% 
  dplyr::filter(X %in% aneuploid_cellline) -> aneuploid_cell_df_flt
#chromosome pos
chr_pos_list = colnames(aneuploid_cell_df_flt)[-1]

# check for the dosage compensation 
dosage_comp = function(x){
  
  # shape of data 
  return_df = lapply(1:dim(x)[1], function(a){
    # gene_name 
    print(a)
    
    gene_name = unlist(x[a,1], use.names = F)
    #chrom_pos
    chr_pos =unlist(x[a,2], use.names = F)
    
    if(chr_pos %in% chr_pos_list){
      # get the corresponding cell lines
      df_flt_anno = aneuploid_cell_df_flt[,c(1,which(colnames(aneuploid_cell_df_flt) == chr_pos))]
      
      # get the loss cells
      loss_cell_name = df_flt_anno[df_flt_anno[,2] == -1,]$X 
      
      # get the neutral cells
      neutral_cell_name = df_flt_anno[df_flt_anno[,2] == 0,]$X 
      
      # get expression
      exp_flt = CCLE_exp[,c(1,which(colnames(CCLE_exp) == gene_name))]
      #loss exp median
      exp_loss = exp_flt[exp_flt[,1] %in% loss_cell_name,][,2]
      # get the median exp
      median_loss_exp = median(exp_loss[!is.na(exp_loss)])
      # neutral exp median
      
      exp_neutral = exp_flt[exp_flt[,1]%in% neutral_cell_name,][,2]
      
      median_neutral_exp = median(exp_neutral[!is.na(exp_neutral)])
      
      result_return = data.frame(gene_name, median_loss_exp, median_neutral_exp)
      
      return(result_return)
      
    }else if(chr_pos == "chrY"){
      # get expression
      exp_flt = CCLE_exp[,c(1,which(colnames(CCLE_exp) == gene_name))]
      #loss exp median
      exp_loss = exp_flt[exp_flt[,1]%in% Y_loss_sample,][,2]
      # get the median exp
      median_loss_exp = median(exp_loss[!is.na(exp_loss)])
      # neutral exp median
      
      exp_neutral = exp_flt[exp_flt[,1]%in% Y_normal,][,2]
      
      median_neutral_exp = median(exp_neutral[!is.na(exp_neutral)])
      
      
    }else{
      
      median_loss_exp = NA
      median_neutral_exp= NA
      result_return = data.frame(gene_name, median_loss_exp, median_neutral_exp)
      
      return(result_return)
    }
    
    
  })
  
  
  return_df_cmb = do.call(rbind, return_df)
  
  return(return_df_cmb)
  
  
}

# sig_gene_exp
sig_gene_exp_df = dosage_comp(sig_loc_gene_df)
# non_sig gene exp
non_sig_gene_exp_df = dosage_comp(non_sig_query_list_flt)

# expression combination
exp_cmb_loc_genes = rbind(sig_gene_exp_df, non_sig_gene_exp_df)

exp_cmb_loc_genes %>% 
  filter(!is.na(median_loss_exp), 
         !is.na(median_neutral_exp)) %>% 
  mutate(norm_exp_para = median_loss_exp - median_neutral_exp,
         norm_exp_control = median_neutral_exp- median_neutral_exp) %>% 
  mutate(condition = ifelse(gene_name %in% sig_gene_exp_df$gene_name, "sig","non_sig")) %>% 
  mutate(condition = factor(condition, levels = c("sig", "non_sig"))) %>% 
  gather(norm_exp_para : norm_exp_control, key= "status", value = "value")-> exp_df_normalised
exp_df_normalised %>% 
  ggplot(aes(x = status, y = value))+
  geom_boxplot(
    aes(color = status, fill = status),
    width = 0.3,
    outlier.shape = NA, outliers = F,
    alpha = 0.4
  ) +
  geom_rangeframe() +
  theme_tufte() +
  scale_fill_manual(values = c( "#7852A9",colors[1])) +
  scale_color_manual(values = c("#7852A9",colors[1] )) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  ylab("Normalised exp")+
  ylim(c(-1,1)) +
  facet_wrap(~condition)-> sig_paralog_compen_plot


# test
# p-value / (2.275e-06)
wilcox.test(exp_df_normalised[exp_df_normalised$condition == "sig" &exp_df_normalised$status == "norm_exp_para"  ,]$value, 
            exp_df_normalised[exp_df_normalised$condition == "non_sig" & exp_df_normalised$status == "norm_exp_para" ,]$value)
#sig median
sig_paralog_median <- median(
  exp_df_normalised$value[
    exp_df_normalised$condition == "sig" &
      exp_df_normalised$status == "norm_exp_para"
  ],
  na.rm = TRUE
)


#non-sig median
nonsig_paralog_median <- median(
  exp_df_normalised$value[
    exp_df_normalised$condition == "non_sig" &
      exp_df_normalised$status == "norm_exp_para"
  ],
  na.rm = TRUE
)
print(paste("sig_paralog median:", sig_paralog_median))
print(paste("nonsig_paralog median:", nonsig_paralog_median))
# plot
ggarrange(sig_paralog_compen_plot, ncol = 2)



############################################
#################Figure 3H #################
############################################
# Project: TCGA paralog gene expression in chromosome loss and neutral cells 

# load the package
library(tidyverse)
library(ggthemes)
library(ggpubr)

# load the data
# Get the TPM data
TCGA_TPM = readRDS("/Volumes/My Passport for Mac/TCGA files/TCGA_TPM_summary_cmb.rds")

#TCGA sample info
TCGA_sample_info = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Aneuploidy patient datasets/mmc2.xlsx",skip = 1)

# get the sig and non-sig paralog genes 
# significant paralog
sig_paralog_df = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/Paralog dependency/Jonathan_analysis/protein_localisation_dataset/sig_37_paralog.xlsx")

# non-significant paralog 
non_sig_paralog_df = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/Paralog dependency/Jonathan_analysis/protein_localisation_dataset/non_sig_paralog.xlsx")

# gene the significant paralog genes on the loss chromosome
# significant genes
sig_paralog_genes = unique(sig_paralog_df$paralog_gene)

# get the non-significant paralog genes
non_sig_paralog_genes = unique(c(non_sig_paralog_df$para_gene_1, non_sig_paralog_df$para_gene_2))
non_sig_paralog_genes_flt = non_sig_paralog_genes[-(non_sig_paralog_genes %in% sig_paralog_genes)]


# get the gene information 
gene_info = readRDS("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/RNAseq/GSE293644_3p_8p/All_gene_info.RDS")
# for the significant query tbl
gene_info %>% 
  dplyr::filter(chromosome_name %in% c(1:22, "X", "Y")) %>% 
  dplyr::filter(hgnc_symbol %in% sig_paralog_genes) %>% 
  dplyr::filter(chromosome_name != "Y") %>% 
  dplyr::mutate(chr_pos = paste(chromosome_name, str_sub(band,1,1), sep = "")) %>% 
  dplyr::select(hgnc_symbol, chr_pos) -> sig_paralog_query_tbl

# for the non-significant query tbl
gene_info %>% 
  dplyr::filter(chromosome_name %in% c(1:22, "X", "Y")) %>% 
  dplyr::filter(hgnc_symbol %in% non_sig_paralog_genes_flt) %>% 
  dplyr::filter(chromosome_name != "Y") %>% 
  dplyr::mutate(chr_pos = paste(chromosome_name, str_sub(band,1,1), sep = "")) %>% 
  dplyr::select(hgnc_symbol, chr_pos) -> non_sig_paralog_query_tbl


#### get the aneuploid sample ####
TCGA_sample_info %>% 
  dplyr::filter(`AneuploidyScore(AS)` > 7) -> TCGA_aneuploid_sample

# rename the sample colnames
trim_TCGA_colname = str_sub(colnames(TCGA_TPM), 1,15)
colnames(TCGA_TPM) = trim_TCGA_colname

######### get a function to compare the relative gene expression ###########
REL_EXP_TCGA = function(tbl){
  # get the query list into the sample
  # the table information
  tbl_shape = dim(tbl)
  
  # gene list
  gene_list = c()
  rel_loss_exp_list = c()
  rel_neutral_exp_list = c()
  for (idx in c(1:tbl_shape[1])) {
    
    # get the gene info
    df = tbl[idx,]
    
    # gene
    gene = df$hgnc_symbol
    # chromosome position 
    chr_pos = df$chr_pos
    
    gene_list[idx] = gene
    if(chr_pos %in% colnames(TCGA_aneuploid_sample)){
      
      # get the sample which chromosome loss 
      x <- TCGA_aneuploid_sample[[chr_pos]]
      
      chr_loss_tbl = TCGA_aneuploid_sample[
        !is.na(x) & x == -1,
        ,
        drop = FALSE
      ]
      
      # get loss sample
      chr_loss_sample = chr_loss_tbl$Sample
      
      # get the sample with chromosome neutral
      
      chr_neutral_tbl = TCGA_aneuploid_sample[
        !is.na(x) & x == 0,
        ,
        drop = FALSE
      ]
      
      # get neutral sample
      chr_neutral_sample = chr_neutral_tbl$Sample
      
      
      # look at the expression and calculate the median expression 
      # loss exp 
      loss_exp = TCGA_TPM[rownames(TCGA_TPM) ==gene, colnames(TCGA_TPM) %in% chr_loss_sample]
      median_loss_exp = median(loss_exp)
      # neutral exp
      neutral_exp = TCGA_TPM[rownames(TCGA_TPM) ==gene, colnames(TCGA_TPM) %in% chr_neutral_sample]
      median_neutral_exp = median(neutral_exp)
      
      # get the relavtive_exp
      # relative loss
      relative_loss = log2(median_loss_exp + 0.5) - log2(median_neutral_exp + 0.5)
      
      # relative neutral
      relative_neutral = log2(median_neutral_exp + 0.5) - log2(median_neutral_exp + 0.5)
      
      
    }else{
      relative_loss = NA
      relative_neutral = NA
    }
    
    # assign the expression
    rel_loss_exp_list[idx] = relative_loss
    rel_neutral_exp_list[idx] = relative_neutral
  }
  
  # generate the return file 
  df_return = data.frame(gene_list, rel_loss_exp_list, rel_neutral_exp_list)
  colnames(df_return) = c("gene", "loss_exp", "neutral_exp")
  
  return(df_return)
}

# sig exp tbl
sig_exp_tbl = REL_EXP_TCGA(sig_paralog_query_tbl)
sig_exp_tbl$condition = "sig"

# non-sig exp tbl
non_sig_exp_tbl = REL_EXP_TCGA(non_sig_paralog_query_tbl)
non_sig_exp_tbl$condition = "non-sig"

# combine the 2 files 
cmb_file = rbind(sig_exp_tbl, non_sig_exp_tbl)

#plot
cmb_file %>%
  dplyr::mutate(condition = factor(condition,
                                   levels = c("sig", "non-sig"))) %>% 
  gather(loss_exp:neutral_exp, 
         key = status,
         value = value) %>% 
  dplyr::mutate(status = factor(status, 
                                levels = c("neutral_exp", "loss_exp"))) -> rel_exp_mdf
rel_exp_mdf %>% 
  ggplot(aes(x = status, y = value))+
  geom_boxplot(aes(color = status, 
                   fill = status), 
               width = 0.3,
               outlier.shape = NA, 
               outliers = F,
               alpha = 0.4) +
  #geom_rangeframe() +
  theme_tufte()+
  scale_color_manual(values = c( "#7852A9",
                                 "#EF8536"))+
  scale_fill_manual(values = c( "#7852A9", 
                                "#EF8536"))+
  labs(x = NULL, y = "Gene expression") +
  facet_wrap(~condition)  -> p

# statistic test
#p-value = 0.004633
wilcox.test(rel_exp_mdf[rel_exp_mdf$condition == "sig" &rel_exp_mdf$status == "loss_exp"  ,]$value, 
            rel_exp_mdf[rel_exp_mdf$condition == "non-sig" & rel_exp_mdf$status == "loss_exp" ,]$value)


ggarrange(p, ncol = 2)


