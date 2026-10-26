##############################################
################## Figure S3 #################
##############################################
# project: RNAi analysis for the paralog dependency


# library
library(tidyverse)
library(MetBrewer)
library(ggthemes)
library(ggpubr)


## load the data

#RNAi dataset
RNAi_df = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /Depmap_RNAi/D2_combined_gene_dep_scores.csv")

# sample information
RNAi_sample_info = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /Depmap_RNAi/sample_info.csv")

# get the CCLE sample info
CCLE_sample_info = read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /Depmap_RNAi/CCLE_sample_info_file_2012-10-18.txt")

# get the Depmap sample info
Depmapsample_info = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /Depmap_RNAi/Model.csv")

#get the colname of RNAi
RNAi_colname = colnames(RNAi_df)

head(RNAi_colname)
# to get the cell line name and its corresponding Depmap id
cell_line_names = c()
for (i in c(2:length(RNAi_colname))) {
  # get the name 
  if(i <= 15){
    name = RNAi_colname[i]
    
    # split the name without "X"
    trimmed_name = str_sub(name, start = 2, end = nchar(name))
    
    cell_line_names = c(cell_line_names,trimmed_name )
  }else{
    trimmed_name = RNAi_colname[i]
    
    cell_line_names = c(cell_line_names,trimmed_name )
    
  }
  
}

# get the cell name to the depmap id
Depmapsample_info_ref = Depmapsample_info[,c("ModelID", "CCLEName")]

cell_line_name_df = as.data.frame(cell_line_names)
colnames(cell_line_name_df) = "CCLEName"
cell_line_name_df %>% 
  left_join(Depmapsample_info_ref, by = "CCLEName") -> join_tbl_CCLE_dep_name

# try to get rid of cell line without the info from Depmap
rm_col_index = as.integer(which(is.na(join_tbl_CCLE_dep_name$ModelID))) + 1

# remove the NA name
RNAi_df_flt = RNAi_df[,-rm_col_index]

# rename the colname into Depmap
colnames(RNAi_df_flt)[2:dim(RNAi_df_flt)[2]] = join_tbl_CCLE_dep_name$ModelID[!is.na(join_tbl_CCLE_dep_name$ModelID)]

# get the gene name and reform 
gene_name = RNAi_df_flt$X
gene_list = c()
for (i in c(1:length(gene_name))) {
  name = gene_name[i]
  
  gene = str_split(name, pattern = " ")[[1]][1]
  
  gene_list = c(gene_list, gene)
}

# rename the X column to the gene name
RNAi_df_flt$X = gene_list


# load the data from aneuploidy score

# dataset
# CCLE aneuploid status annotation
aneuploid_dep = read.csv("/Volumes/My Passport for Mac/Postdoc_sheltzer/Depmp data /aneuploid depmap/aneuploidy_scores.csv")
# CCLE cell line with each chromsome arm gain or loss info
df_aneu_score = read.csv("/Volumes/My Passport for Mac/Postdoc_sheltzer/Depmp data /aneuploid depmap/arm_call_scores.csv")
#cell line annotation
cell_line_id = read.delim("/Volumes/My Passport for Mac/Postdoc_sheltzer/Depmp data /aneuploid depmap/Cell_lines_annotations_20181226.txt")

# color set
colors = c("#EF8536","#3A76AF")


# correlation dependency
## analysis here: correlation and dep and exp
# retrieve gene expression data
CCLE_exp = read.csv("/Volumes/My Passport for Mac/Postdoc_sheltzer/Depmp data /OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")

# change the name of CCL2_exp data
colnames(CCLE_exp)[1] = "DepMap_ID"
#dim(df_gene_dep)

# gene the gene name 
#dim(df_gene_dep)
genelist_exp = sub("\\..*","",colnames(CCLE_exp)[2:dim(CCLE_exp)[2]])

# replace gene to colnames
colnames(CCLE_exp)[2:dim(CCLE_exp)[2]] = genelist_exp


##################### 8p loss sample ##################

# definition: 
# Aneuploidy score  median = 25: high aneuploidy
# Near-euploid; median=3
# to avoid the bias from other chromosome and get accuracy dependecy on 7p gain
# we have to set up an aneuploidy score cut-off

#get the 8p loss cells
df_aneu_score %>% 
  filter(X8p == -1) %>% 
  dplyr::select(X) %>%
  unlist(use.names = F) -> chr_8p_loss

df_aneu_score %>% 
  filter(X8p == 0) %>% 
  dplyr::select(X) %>% 
  unlist(use.names = F) -> chr_8p_disomy


# Primary setting
# based on the aneuploidy score
# 8p monosomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_8p_loss) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 27) %>%
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_8p_loss_pri_flt

# 1q disomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_8p_disomy) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 27) %>% 
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_8p_diso_pri_flt  


# compare the RNAi dependency

RNAi_df_flt %>% 
  filter(X == "PPP2CA") -> PPP2CA_dep_RNAi

PPP2CA_dep_RNAi_t = as.data.frame(t(PPP2CA_dep_RNAi))

# rename the dep 
colnames(PPP2CA_dep_RNAi_t) = c("PPP2CA_dep")

PPP2CA_dep_RNAi_t$DepID = rownames(PPP2CA_dep_RNAi_t)
PPP2CA_dep_RNAi_t_edt = PPP2CA_dep_RNAi_t[-c(1),]

# assign the 8p loss and no loss 
PPP2CA_dep_RNAi_t_edt %>% 
  dplyr::mutate(PPP2CA_dep = as.numeric(PPP2CA_dep)) %>% 
  filter(DepID %in% c(chr_8p_loss_pri_flt, chr_8p_diso_pri_flt)) %>% 
  dplyr::mutate(condition = ifelse(DepID %in%chr_8p_loss_pri_flt ,"8p_loss", "8p_non_loss")) -> PPP2CA_dep

# plot the dependency 
PPP2CA_dep %>% 
  ggplot(aes(x = condition, y = PPP2CA_dep)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  theme(legend.position = "none") +
  labs(y = "PPP2CA_dep (RNAi)")

# t test 
# p value = p-value = 3.987e-05
t.test(PPP2CA_dep[PPP2CA_dep$condition == "8p_loss", ]$PPP2CA_dep, PPP2CA_dep[PPP2CA_dep$condition == "8p_non_loss", ]$PPP2CA_dep)

# test for the correlation between 

# dependency and exp correlation analysis
PPP2CA_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  PPP2CA_dep %>% 
    rename(DepMap_ID = DepID) %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,1],comb_tmp[,4])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})


# combine the list
PPP2CA_exp_corr_cmb = do.call( rbind,PPP2CA_exp_corr)
colnames(PPP2CA_exp_corr_cmb) = c("gene", "p_value", "corr")


library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
PPP2CA_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(PPP2CA_exp_corr_cmb,gene%in% c("PPP2CB")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    #label.padding = 0.25,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
  ) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[2],"grey",colors[1])) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         legend.position = "none")


# correlation plot
# based on the expression and correlation of PPP2CB

PPP2CB_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "PPP2CB"))]
PPP2CA_dep %>% 
  rename(DepMap_ID = DepID) %>% 
  left_join(PPP2CB_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(PPP2CA_dep)) %>% 
  filter(!is.na(PPP2CB)) -> cmb_df
cmb_df$pc = predict(prcomp(~PPP2CA_dep+PPP2CB, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = PPP2CB, y = PPP2CA_dep)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = -1)) +
  theme(legend.position = "none") +
  labs( y = "PPP2CA_dep (RNAi)", x = "PPP2CB exp (TPM)") +
  geom_text( x = 7, y =-1.58 , label = "R = 0.365549")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value = 8.946e-12", x = 6.8, y =-1.7) 

# correlation test
#p-value = 8.946e-12
cor_result = cor.test(cmb_df$PPP2CA_dep, cmb_df$PPP2CB)



################### Y chr loss ########################

# gene copy number
CCLE_gene_cnv = read.csv("/Volumes/My Passport for Mac/Postdoc_sheltzer/Depmp data /PortalOmicsCNGeneLog2.csv")

# cell and patient info
cell_line_id = read.delim("/Volumes/My Passport for Mac/Postdoc_sheltzer/Depmp data /aneuploid depmap/Cell_lines_annotations_20181226.txt")

# rename the gene name on CCLE_gene_CNV
# change the name of data
colnames(CCLE_gene_cnv)[1] = "DepMap_ID"
#dim(df_gene_dep)

# gene the gene name 
dim(CCLE_gene_cnv)
genelist_Dep = sub("\\..*","",colnames(CCLE_gene_cnv)[2:dim(CCLE_gene_cnv)[2]])

# replace gene to colnames
colnames(CCLE_gene_cnv)[2:dim(CCLE_gene_cnv)[2]] = genelist_Dep
head(CCLE_gene_cnv[,1:6])

library(biomaRt)
# retrieve Y chromosome genes
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
results <- getBM(attributes = c("chromosome_name", "entrezgene_id", "hgnc_symbol"),
                 filters = "chromosome_name", values = "Y", mart = mart)

Y_chr_genename = results$hgnc_symbol[results$hgnc_symbol != ""]

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
Y_normal = Y_cnv_male$DepMap_ID[!(Y_cnv_male$DepMap_ID %in% Y_loss_sample)]



# The gene we found in Y chr loss dependency is EIF1AX and 
RNAi_df_flt %>% 
  .[1151,]  -> EIF1AX_dep_RNAi

EIF1AX_dep_RNAi_t = as.data.frame(t(EIF1AX_dep_RNAi))

# rename the dep 
colnames(EIF1AX_dep_RNAi_t) = c("EIF1AX_dep")

EIF1AX_dep_RNAi_t$DepID = rownames(EIF1AX_dep_RNAi_t)
EIF1AX_dep_RNAi_t_edt = EIF1AX_dep_RNAi_t[-c(1),]

# assign the Y loss and no loss 
EIF1AX_dep_RNAi_t_edt %>% 
  dplyr::mutate(EIF1AX_dep = as.numeric(EIF1AX_dep)) %>% 
  filter(DepID %in% c(Y_loss_sample, Y_normal)) %>% 
  dplyr::mutate(condition = ifelse(DepID %in%Y_loss_sample ,"Y_loss", "Y_non_loss")) -> EIF1AX_dep

# plot the dependency 
EIF1AX_dep %>% 
  ggplot(aes(x = condition, y = EIF1AX_dep)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  theme(legend.position = "none") +
  labs(y = "EIF1AX_dep (RNAi)")

# t test 
#  p-value = 0.002505
# only a few cell line in the RNAi
t.test(EIF1AX_dep[EIF1AX_dep$condition == "Y_loss", ]$EIF1AX_dep, EIF1AX_dep[EIF1AX_dep$condition == "Y_non_loss", ]$EIF1AX_dep)

# test for the correlation between 

# dependency and exp correlation analysis
EIF1AX_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  EIF1AX_dep %>% 
    rename(DepMap_ID = DepID) %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,1],comb_tmp[,4])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})


# combine the list
EIF1AX_exp_corr_cmb = do.call( rbind,EIF1AX_exp_corr)
colnames(EIF1AX_exp_corr_cmb) = c("gene", "p_value", "corr")


library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
EIF1AX_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(EIF1AX_exp_corr_cmb,gene%in% c("EIF1AY")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    #label.padding = 0.25,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
  ) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[2],"grey",colors[1])) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         legend.position = "none")


# correlation plot
# based on the expression and correlation of PPP2CB

EIF1AY_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "EIF1AY"))]
EIF1AX_dep %>% 
  rename(DepMap_ID = DepID) %>% 
  left_join(EIF1AY_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(EIF1AX_dep)) %>% 
  filter(!is.na(EIF1AY)) -> cmb_df
cmb_df$pc = predict(prcomp(~EIF1AX_dep+EIF1AY, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = EIF1AY, y = EIF1AX_dep)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = -1)) +
  theme(legend.position = "none") +
  labs( y = "EIF1AX_dep (RNAi)", x = "EIF1AY exp (TPM)") +
  geom_text( x = 2, y = -0.2 , label = "R = 0.4051632 ")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value = 6.113e-10", x = 2, y =-0.25) 

# correlation test
#p-value = 6.113e-10
cor_result = cor.test(cmb_df$EIF1AX_dep, cmb_df$EIF1AY)


# RPS4X ###

# The gene we found in Y chr loss dependency is RPS4X and 
RNAi_df_flt %>% 
  filter(X == "RPS4X")  -> RPS4X_dep_RNAi

RPS4X_dep_RNAi_t = as.data.frame(t(RPS4X_dep_RNAi))

# rename the dep 
colnames(RPS4X_dep_RNAi_t) = c("RPS4X_dep")

RPS4X_dep_RNAi_t$DepID = rownames(RPS4X_dep_RNAi_t)
RPS4X_dep_RNAi_t_edt = RPS4X_dep_RNAi_t[-c(1),]

# assign the Y loss and no loss 
RPS4X_dep_RNAi_t_edt %>% 
  dplyr::mutate(RPS4X_dep = as.numeric(RPS4X_dep)) %>% 
  filter(DepID %in% c(Y_loss_sample, Y_normal)) %>% 
  dplyr::mutate(condition = ifelse(DepID %in%Y_loss_sample ,"Y_loss", "Y_non_loss")) -> RPS4X_dep

# plot the dependency 
RPS4X_dep %>% 
  ggplot(aes(x = condition, y = RPS4X_dep)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  theme(legend.position = "none") +
  labs(y = "RPS4X_dep (RNAi)")

# t test 
#  p-value = 1.43e-05
# only a few cell line in the RNAi
t.test(RPS4X_dep[RPS4X_dep$condition == "Y_loss", ]$RPS4X_dep, RPS4X_dep[RPS4X_dep$condition == "Y_non_loss", ]$RPS4X_dep)

# test for the correlation between 

# dependency and exp correlation analysis
RPS4X_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  RPS4X_dep %>% 
    rename(DepMap_ID = DepID) %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,1],comb_tmp[,4])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})


# combine the list
RPS4X_exp_corr_cmb = do.call( rbind,RPS4X_exp_corr)
colnames(RPS4X_exp_corr_cmb) = c("gene", "p_value", "corr")


library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
RPS4X_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(RPS4X_exp_corr_cmb,gene%in% c("RPS4Y1")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    #label.padding = 0.25,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
  ) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[2],"grey",colors[1])) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         legend.position = "none")


# correlation plot
# based on the expression and correlation of PPP2CB

RPS4Y1_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "RPS4Y1"))]
RPS4X_dep %>% 
  rename(DepMap_ID = DepID) %>% 
  left_join(RPS4Y1_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(RPS4X_dep)) %>% 
  filter(!is.na(RPS4Y1)) -> cmb_df
cmb_df$pc = predict(prcomp(~RPS4X_dep+RPS4Y1, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = RPS4Y1, y = RPS4X_dep)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = -1)) +
  theme(legend.position = "none") +
  labs( y = "RPS4X_dep (RNAi)", x = "RPS4Y1 exp (TPM)") +
  geom_text( x = 2, y = -0.2 , label = "R = 0.5176328 ")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value = 6.113e-10", x = 2, y =-0.3) 

# correlation test
#p-value = 6.113e-10
cor_result = cor.test(cmb_df$RPS4X_dep, cmb_df$RPS4Y1)



################### 9p loss #####################
#get the 9p loss cells
df_aneu_score %>% 
  filter(X9p == -1) %>% 
  dplyr::select(X) %>%
  unlist(use.names = F) -> chr_9p_loss

df_aneu_score %>% 
  filter(X9p == 0) %>% 
  dplyr::select(X) %>% 
  unlist(use.names = F) -> chr_9p_disomy


# Primary setting
# based on the aneuploidy score
# 9p monosomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_9p_loss) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 26) %>%
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_9p_loss_pri_flt

# 9p disomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_9p_disomy) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 26) %>% 
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_9p_diso_pri_flt  



RNAi_df_flt %>% 
  filter(X == "SMARCA4") -> SMARCA4_dep_RNAi

SMARCA4_dep_RNAi_t = as.data.frame(t(SMARCA4_dep_RNAi))

# rename the dep 
colnames(SMARCA4_dep_RNAi_t) = c("SMARCA4_dep")

SMARCA4_dep_RNAi_t$DepID = rownames(SMARCA4_dep_RNAi_t)
SMARCA4_dep_RNAi_t_edt = SMARCA4_dep_RNAi_t[-c(1),]

#View(SMARCA4_dep_RNAi_t_edt)
# assign the 9p loss and no loss 
SMARCA4_dep_RNAi_t_edt %>% 
  dplyr::mutate(SMARCA4_dep = as.numeric(SMARCA4_dep)) %>% 
  filter(DepID %in% c(chr_9p_loss_pri_flt, chr_9p_diso_pri_flt)) %>% 
  dplyr::mutate(condition = ifelse(DepID %in%chr_9p_loss_pri_flt ,"9p_loss", "9p_non_loss")) -> SMARCA4_dep

# plot the dependency 
SMARCA4_dep %>% 
  ggplot(aes(x = condition, y = SMARCA4_dep)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  theme(legend.position = "none") +
  labs(y = "SMARCA4_dep (RNAi)")


# only a few cell line in the RNAi
t.test(SMARCA4_dep[SMARCA4_dep$condition == "9p_loss",]$SMARCA4_dep, SMARCA4_dep[SMARCA4_dep$condition == "9p_non_loss",]$SMARCA4_dep)

# test for the correlation between 

# dependency and exp correlation analysis
SMARCA4_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  SMARCA4_dep %>% 
    rename(DepMap_ID = DepID) %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,1],comb_tmp[,4])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})


# combine the list
SMARCA4_exp_corr_cmb = do.call( rbind,SMARCA4_exp_corr)
colnames(SMARCA4_exp_corr_cmb) = c("gene", "p_value", "corr")


library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
SMARCA4_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(SMARCA4_exp_corr_cmb,gene%in% c("SMARCA2")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    #label.padding = 0.25,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
  ) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[2],"grey",colors[1])) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         legend.position = "none")



SMARCA2_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "SMARCA2"))]
SMARCA4_dep %>% 
  rename(DepMap_ID = DepID) %>% 
  left_join(SMARCA2_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(SMARCA4_dep)) %>% 
  filter(!is.na(SMARCA2)) -> cmb_df
cmb_df$pc = predict(prcomp(~SMARCA4_dep+SMARCA2, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = SMARCA2, y = SMARCA4_dep)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = -1)) +
  theme(legend.position = "none") +
  labs( y = "SMARCA4_dep (RNAi)", x = "SMARCA2 exp (TPM)") +
  geom_text( x = 1.4, y = 0.58 , label = "R = 0.03977871  ")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value = 0.4714", x = 1.4, y =0.5) 

# correlation test
#p-value = 0.4714
cor_result = cor.test(cmb_df$SMARCA4_dep, cmb_df$SMARCA2)



################## 6p loss ##################
#get the 6p loss cells
df_aneu_score %>% 
  filter(X6p == -1) %>% 
  dplyr::select(X) %>%
  unlist(use.names = F) -> chr_6p_loss

df_aneu_score %>% 
  filter(X6p == 0) %>% 
  dplyr::select(X) %>% 
  unlist(use.names = F) -> chr_6p_disomy


# Primary setting
# based on the aneuploidy score
# 6p monosomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_6p_loss) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 26) %>%
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_6p_loss_pri_flt

# 6p disomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_6p_disomy) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 26) %>% 
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_6p_diso_pri_flt  



# RNAi dependency
RNAi_df_flt %>% 
  filter(X == "TUBB4B") -> TUBB4B_dep_RNAi

TUBB4B_dep_RNAi_t = as.data.frame(t(TUBB4B_dep_RNAi))

# rename the dep 
colnames(TUBB4B_dep_RNAi_t) = c("TUBB4B_dep")

TUBB4B_dep_RNAi_t$DepID = rownames(TUBB4B_dep_RNAi_t)
TUBB4B_dep_RNAi_t_edt = TUBB4B_dep_RNAi_t[-c(1),]


# assign the 6p loss and no loss 
TUBB4B_dep_RNAi_t_edt %>% 
  dplyr::mutate(TUBB4B_dep = as.numeric(TUBB4B_dep)) %>% 
  filter(DepID %in% c(chr_6p_loss_pri_flt, chr_6p_diso_pri_flt)) %>% 
  dplyr::mutate(condition = ifelse(DepID %in%chr_6p_loss_pri_flt ,"6p_loss", "6p_non_loss")) -> TUBB4B_dep

# plot the dependency 
TUBB4B_dep %>% 
  ggplot(aes(x = condition, y = TUBB4B_dep)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  theme(legend.position = "none") +
  labs(y = "TUBB4B_dep (RNAi)")

# p value = 0.01573
# only a few cell line in the RNAi
t.test(TUBB4B_dep[TUBB4B_dep$condition == "6p_loss",]$TUBB4B_dep, TUBB4B_dep[TUBB4B_dep$condition == "6p_non_loss",]$TUBB4B_dep)

# test for the correlation between 

# dependency and exp correlation analysis
TUBB4B_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  TUBB4B_dep %>% 
    rename(DepMap_ID = DepID) %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,1],comb_tmp[,4])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})


# combine the list
TUBB4B_exp_corr_cmb = do.call( rbind,TUBB4B_exp_corr)
colnames(TUBB4B_exp_corr_cmb) = c("gene", "p_value", "corr")


library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
TUBB4B_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(TUBB4B_exp_corr_cmb,gene%in% c("TUBB", "TUBB2B", "TUBB4A", "TUBB2A")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    #label.padding = 0.25,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
  ) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[2],"grey",colors[1])) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         legend.position = "none")



TUBB_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "TUBB"))]
TUBB4B_dep %>% 
  rename(DepMap_ID = DepID) %>% 
  left_join(TUBB_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(TUBB4B_dep)) %>% 
  filter(!is.na(TUBB)) -> cmb_df
cmb_df$pc = predict(prcomp(~TUBB4B_dep+TUBB, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = TUBB, y = TUBB4B_dep)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = -1)) +
  theme(legend.position = "none") +
  labs( y = "TUBB4B_dep (RNAi)", x = "TUBB exp (TPM)") +
  geom_text( x = 10, y = -0.5, label = "R = 0.2864524")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value = 6.799e-06", x = 10, y =-0.55) 

# correlation test
#p-value = 2.01e-08
cor_result = cor.test(cmb_df$TUBB4B_dep, cmb_df$TUBB)

################### 11q loss ######################

#get the 11q loss cells
df_aneu_score %>% 
  filter(X11q == -1) %>% 
  dplyr::select(X) %>%
  unlist(use.names = F) -> chr_11q_loss

df_aneu_score %>% 
  filter(X11q == 0) %>% 
  dplyr::select(X) %>% 
  unlist(use.names = F) -> chr_11q_disomy


# Primary setting
# based on the aneuploidy score
# 11q monosomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_11q_loss) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 27) %>%
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_11q_loss_pri_flt

# 11q disomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_11q_disomy) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 27) %>% 
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_11q_diso_pri_flt  


# RNAi dependency
RNAi_df_flt %>% 
  filter(X == "PPP2R1A") -> PPP2R1A_dep_RNAi

PPP2R1A_dep_RNAi_t = as.data.frame(t(PPP2R1A_dep_RNAi))

# rename the dep 
colnames(PPP2R1A_dep_RNAi_t) = c("PPP2R1A_dep")

PPP2R1A_dep_RNAi_t$DepID = rownames(PPP2R1A_dep_RNAi_t)
PPP2R1A_dep_RNAi_t_edt = PPP2R1A_dep_RNAi_t[-c(1),]


# assign the 11q loss and no loss 
PPP2R1A_dep_RNAi_t_edt %>% 
  dplyr::mutate(PPP2R1A_dep = as.numeric(PPP2R1A_dep)) %>% 
  filter(DepID %in% c(chr_11q_loss_pri_flt, chr_11q_diso_pri_flt)) %>% 
  dplyr::mutate(condition = ifelse(DepID %in%chr_11q_loss_pri_flt ,"11q_loss", "11q_non_loss")) -> PPP2R1A_dep

# plot the dependency 
PPP2R1A_dep %>% 
  ggplot(aes(x = condition, y = PPP2R1A_dep)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  theme(legend.position = "none") +
  labs(y = "PPP2R1A_dep (RNAi)")

# p value = 2.192e-06
# only a few cell line in the RNAi
t.test(PPP2R1A_dep[PPP2R1A_dep$condition == "11q_loss",]$PPP2R1A_dep, PPP2R1A_dep[PPP2R1A_dep$condition == "11q_non_loss",]$PPP2R1A_dep)

# test for the correlation between 

# dependency and exp correlation analysis
PPP2R1A_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  PPP2R1A_dep %>% 
    rename(DepMap_ID = DepID) %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,1],comb_tmp[,4])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})


# combine the list
PPP2R1A_exp_corr_cmb = do.call( rbind,PPP2R1A_exp_corr)
colnames(PPP2R1A_exp_corr_cmb) = c("gene", "p_value", "corr")


library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
PPP2R1A_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(PPP2R1A_exp_corr_cmb,gene%in% c("PPP2R1B")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    #label.padding = 0.25,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
  ) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[2],"grey",colors[1])) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         legend.position = "none")



PPP2R1B_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "PPP2R1B"))]
PPP2R1A_dep %>% 
  rename(DepMap_ID = DepID) %>% 
  left_join(PPP2R1B_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(PPP2R1A_dep)) %>% 
  filter(!is.na(PPP2R1B)) -> cmb_df
cmb_df$pc = predict(prcomp(~PPP2R1A_dep+PPP2R1B, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = PPP2R1B, y = PPP2R1A_dep)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = -1)) +
  theme(legend.position = "none") +
  labs( y = "PPP2R1A_dep (RNAi)", x = "PPP2R1B exp (TPM)") +
  geom_text( x = 6, y = -2, label = "R = 0.3134336 ")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value =2.007e-08", x = 6 ,y =-2.2) 

# correlation test
#p-value = 2.01e-08
cor_result = cor.test(cmb_df$PPP2R1A_dep, cmb_df$PPP2R1B)


################# 10q ###################
#get the 10q loss cells
df_aneu_score %>% 
  filter(X10q == -1) %>% 
  dplyr::select(X) %>%
  unlist(use.names = F) -> chr_10q_loss

df_aneu_score %>% 
  filter(X10q == 0) %>% 
  dplyr::select(X) %>% 
  unlist(use.names = F) -> chr_10q_disomy


# Primary setting
# based on the aneuploidy score
# 10q monosomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_10q_loss) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 26.5) %>%
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_10q_loss_pri_flt

# 10q disomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_10q_disomy) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 26.5) %>% 
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_10q_diso_pri_flt  

# RNAi dependency
RNAi_df_flt %>% 
  filter(X == "FBXW11") -> FBXW11_dep_RNAi

FBXW11_dep_RNAi_t = as.data.frame(t(FBXW11_dep_RNAi))

# rename the dep 
colnames(FBXW11_dep_RNAi_t) = c("FBXW11_dep")

FBXW11_dep_RNAi_t$DepID = rownames(FBXW11_dep_RNAi_t)
FBXW11_dep_RNAi_t_edt = FBXW11_dep_RNAi_t[-c(1),]


# assign the 10q loss and no loss 
FBXW11_dep_RNAi_t_edt %>% 
  dplyr::mutate(FBXW11_dep = as.numeric(FBXW11_dep)) %>% 
  filter(DepID %in% c(chr_10q_loss_pri_flt, chr_10q_diso_pri_flt)) %>% 
  dplyr::mutate(condition = ifelse(DepID %in%chr_10q_loss_pri_flt ,"10q_loss", "10q_non_loss")) -> FBXW11_dep

# plot the dependency 
FBXW11_dep %>% 
  ggplot(aes(x = condition, y = FBXW11_dep)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  theme(legend.position = "none") +
  labs(y = "FBXW11_dep (RNAi)")

# p value = 0.01586
# only a few cell line in the RNAi
t.test(FBXW11_dep[FBXW11_dep$condition == "10q_loss",]$FBXW11_dep, FBXW11_dep[FBXW11_dep$condition == "10q_non_loss",]$FBXW11_dep)

# test for the correlation between 

# dependency and exp correlation analysis
FBXW11_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  FBXW11_dep %>% 
    rename(DepMap_ID = DepID) %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,1],comb_tmp[,4])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})


# combine the list
FBXW11_exp_corr_cmb = do.call( rbind,FBXW11_exp_corr)
colnames(FBXW11_exp_corr_cmb) = c("gene", "p_value", "corr")


library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
FBXW11_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(FBXW11_exp_corr_cmb,gene%in% c("BTRC")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    #label.padding = 0.25,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
  ) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[2],"grey",colors[1])) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         legend.position = "none")



BTRC_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "BTRC"))]
FBXW11_dep %>% 
  rename(DepMap_ID = DepID) %>% 
  left_join(BTRC_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(FBXW11_dep)) %>% 
  filter(!is.na(BTRC)) -> cmb_df
cmb_df$pc = predict(prcomp(~FBXW11_dep+BTRC, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = BTRC, y = FBXW11_dep)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = 1)) +
  theme(legend.position = "none") +
  labs( y = "FBXW11_dep (RNAi)", x = "BTRC exp (TPM)") +
  geom_text( x = 4.2, y = -0.7, label = "R = 0.1933132 ")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value =0.0002862", x = 4.2 ,y =-0.78) 

# correlation test
#p-value = 0.0002862
cor_result = cor.test(cmb_df$FBXW11_dep, cmb_df$BTRC)



####################### 19p ######################
#get the 19p loss cells
df_aneu_score %>% 
  filter(X19p == -1) %>% 
  dplyr::select(X) %>%
  unlist(use.names = F) -> chr_19p_loss

df_aneu_score %>% 
  filter(X19p == 0) %>% 
  dplyr::select(X) %>% 
  unlist(use.names = F) -> chr_19p_disomy


# Primary setting
# based on the aneuploidy score
# 19p monosomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_19p_loss) %>% 
  filter(Aneuploidy.score > 8 & Aneuploidy.score < 25) %>%
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_19p_loss_pri_flt

# 19p disomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_19p_disomy) %>% 
  filter(Aneuploidy.score > 8 & Aneuploidy.score < 25) %>% 
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_19p_diso_pri_flt  


# RNAi dependency
RNAi_df_flt %>% 
  filter(X == "DDX39B") -> DDX39B_dep_RNAi

DDX39B_dep_RNAi_t = as.data.frame(t(DDX39B_dep_RNAi))

# rename the dep 
colnames(DDX39B_dep_RNAi_t) = c("DDX39B_dep")

DDX39B_dep_RNAi_t$DepID = rownames(DDX39B_dep_RNAi_t)
DDX39B_dep_RNAi_t_edt = DDX39B_dep_RNAi_t[-c(1),]


# assign the 19p loss and no loss 
DDX39B_dep_RNAi_t_edt %>% 
  dplyr::mutate(DDX39B_dep = as.numeric(DDX39B_dep)) %>% 
  filter(DepID %in% c(chr_19p_loss_pri_flt, chr_19p_diso_pri_flt)) %>% 
  dplyr::mutate(condition = ifelse(DepID %in%chr_19p_loss_pri_flt ,"19p_loss", "19p_non_loss")) -> DDX39B_dep

# plot the dependency 
DDX39B_dep %>% 
  ggplot(aes(x = condition, y = DDX39B_dep)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  theme(legend.position = "none") +
  labs(y = "DDX39B_dep (RNAi)")

# p value =  5.519e-06
# only a few cell line in the RNAi
t.test(DDX39B_dep[DDX39B_dep$condition == "19p_loss",]$DDX39B_dep, DDX39B_dep[DDX39B_dep$condition == "19p_non_loss",]$DDX39B_dep)

# test for the correlation between 

# dependency and exp correlation analysis
DDX39B_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  DDX39B_dep %>% 
    rename(DepMap_ID = DepID) %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,1],comb_tmp[,4])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})


# combine the list
DDX39B_exp_corr_cmb = do.call( rbind,DDX39B_exp_corr)
colnames(DDX39B_exp_corr_cmb) = c("gene", "p_value", "corr")


library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
DDX39B_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(DDX39B_exp_corr_cmb,gene%in% c("DDX39A")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    #label.padding = 0.25,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
  ) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[2],"grey",colors[1])) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         legend.position = "none")



DDX39A_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "DDX39A"))]
DDX39B_dep %>% 
  rename(DepMap_ID = DepID) %>% 
  left_join(DDX39A_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(DDX39B_dep)) %>% 
  filter(!is.na(DDX39A)) -> cmb_df
cmb_df$pc = predict(prcomp(~DDX39B_dep+DDX39A, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = DDX39A, y = DDX39B_dep)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = 1)) +
  theme(legend.position = "none") +
  labs( y = "DDX39B_dep (RNAi)", x = "DDX39A exp (TPM)") +
  geom_text( x = 7.5, y = -1.5, label = "R = 0.3654578 ")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value =5.651e-11", x = 7.5 ,y =-1.58) 

# correlation test
#p-value = 5.651e-11
cor_result = cor.test(cmb_df$DDX39B_dep, cmb_df$DDX39A)


####################### 12p loss ########################
#get the 12p loss cells
df_aneu_score %>% 
  filter(X12p == -1) %>% 
  dplyr::select(X) %>%
  unlist(use.names = F) -> chr_12p_loss

df_aneu_score %>% 
  filter(X12p == 0) %>% 
  dplyr::select(X) %>% 
  unlist(use.names = F) -> chr_12p_disomy


# Primary setting
# based on the aneuploidy score
# 12p monosomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_12p_loss) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 26) %>%
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_12p_loss_pri_flt

# 15q disomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_12p_disomy) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 26) %>% 
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_12p_diso_pri_flt  



# RNAi dependency
RNAi_df_flt %>% 
  filter(X == "MAGOH") -> MAGOH_dep_RNAi

MAGOH_dep_RNAi_t = as.data.frame(t(MAGOH_dep_RNAi))

# rename the dep 
colnames(MAGOH_dep_RNAi_t) = c("MAGOH_dep")

MAGOH_dep_RNAi_t$DepID = rownames(MAGOH_dep_RNAi_t)
MAGOH_dep_RNAi_t_edt = MAGOH_dep_RNAi_t[-c(1),]


# assign the 12p loss and no loss 
MAGOH_dep_RNAi_t_edt %>% 
  dplyr::mutate(MAGOH_dep = as.numeric(MAGOH_dep)) %>% 
  filter(DepID %in% c(chr_12p_loss_pri_flt, chr_12p_diso_pri_flt)) %>% 
  dplyr::mutate(condition = ifelse(DepID %in%chr_12p_loss_pri_flt ,"12p_loss", "12p_non_loss")) -> MAGOH_dep

# plot the dependency 
MAGOH_dep %>% 
  ggplot(aes(x = condition, y = MAGOH_dep)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  theme(legend.position = "none") +
  labs(y = "MAGOH_dep (RNAi)")

# p value =  5.124e-07
# only a few cell line in the RNAi
t.test(MAGOH_dep[MAGOH_dep$condition == "12p_loss",]$MAGOH_dep, MAGOH_dep[MAGOH_dep$condition == "12p_non_loss",]$MAGOH_dep)

# test for the correlation between 

# dependency and exp correlation analysis
MAGOH_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  MAGOH_dep %>% 
    rename(DepMap_ID = DepID) %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,1],comb_tmp[,4])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})


# combine the list
MAGOH_exp_corr_cmb = do.call( rbind,MAGOH_exp_corr)
colnames(MAGOH_exp_corr_cmb) = c("gene", "p_value", "corr")


library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
MAGOH_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(MAGOH_exp_corr_cmb,gene%in% c("MAGOHB")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    #label.padding = 0.25,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
  ) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[2],"grey",colors[1])) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         legend.position = "none")



MAGOHB_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "MAGOHB"))]
MAGOH_dep %>% 
  rename(DepMap_ID = DepID) %>% 
  left_join(MAGOHB_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(MAGOH_dep)) %>% 
  filter(!is.na(MAGOHB)) -> cmb_df
cmb_df$pc = predict(prcomp(~MAGOH_dep+MAGOHB, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = MAGOHB, y = MAGOH_dep)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = 1)) +
  theme(legend.position = "none") +
  labs( y = "MAGOH_dep (RNAi)", x = "MAGOHB exp (TPM)") +
  geom_text( x = 5.8, y = -0.8, label = "R = 0.3675493")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value =4.07e-08", x = 5.8 ,y =-0.9) 

# correlation test
#p-value = 4.07e-08
cor_result = cor.test(cmb_df$MAGOH_dep, cmb_df$MAGOHB)


##################### 3p loss #######################
#get the 3p loss cells
df_aneu_score %>% 
  filter(X3p == -1) %>% 
  dplyr::select(X) %>%
  unlist(use.names = F) -> chr_3p_loss

df_aneu_score %>% 
  filter(X3p == 0) %>% 
  dplyr::select(X) %>% 
  unlist(use.names = F) -> chr_3p_disomy



# Primary setting
# based on the aneuploidy score
# 3p monosomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_3p_loss) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 26) %>%
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_3p_loss_pri_flt

# 3p disomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_3p_disomy) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 26) %>% 
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_3p_diso_pri_flt 

# RNAi dependency
RNAi_df_flt %>% 
  filter(X == "ARF1") -> ARF1_dep_RNAi

ARF1_dep_RNAi_t = as.data.frame(t(ARF1_dep_RNAi))

# rename the dep 
colnames(ARF1_dep_RNAi_t) = c("ARF1_dep")

ARF1_dep_RNAi_t$DepID = rownames(ARF1_dep_RNAi_t)
ARF1_dep_RNAi_t_edt = ARF1_dep_RNAi_t[-c(1),]


# assign the 3pp loss and no loss 
ARF1_dep_RNAi_t_edt %>% 
  dplyr::mutate(ARF1_dep = as.numeric(ARF1_dep)) %>% 
  filter(DepID %in% c(chr_3p_loss_pri_flt, chr_3p_diso_pri_flt)) %>% 
  dplyr::mutate(condition = ifelse(DepID %in%chr_3p_loss_pri_flt ,"3p_loss", "3p_non_loss")) -> ARF1_dep

# plot the dependency 
ARF1_dep %>% 
  ggplot(aes(x = condition, y = ARF1_dep)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  theme(legend.position = "none") +
  labs(y = "ARF1_dep (RNAi)")

# p value =  0.05163
# only a few cell line in the RNAi
t.test(ARF1_dep[ARF1_dep$condition == "3p_loss",]$ARF1_dep, ARF1_dep[ARF1_dep$condition == "3p_non_loss",]$ARF1_dep)

# test for the correlation between 

# dependency and exp correlation analysis
ARF1_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  ARF1_dep %>% 
    rename(DepMap_ID = DepID) %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,1],comb_tmp[,4])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})


# combine the list
ARF1_exp_corr_cmb = do.call( rbind,ARF1_exp_corr)
colnames(ARF1_exp_corr_cmb) = c("gene", "p_value", "corr")


library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
ARF1_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(ARF1_exp_corr_cmb,gene%in% c("ARF4")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    #label.padding = 0.25,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
  ) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[2],"grey",colors[1])) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         legend.position = "none")



ARF4_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "ARF4"))]
ARF1_dep %>% 
  rename(DepMap_ID = DepID) %>% 
  left_join(ARF4_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(ARF1_dep)) %>% 
  filter(!is.na(ARF4)) -> cmb_df
cmb_df$pc = predict(prcomp(~ARF1_dep+ARF4, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = ARF4, y = ARF1_dep)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = 1)) +
  theme(legend.position = "none") +
  labs( y = "ARF1_dep (RNAi)", x = "ARF4 exp (TPM)") +
  geom_text( x = 5.8, y = -0.8, label = "R = 0.03881449 ")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value =0.4769", x = 5.8 ,y =-0.9) 

# correlation test
#p-value = 0.4769
cor_result = cor.test(cmb_df$ARF1_dep, cmb_df$ARF4)



################# chr22q ##################
#get the 22q loss cells
df_aneu_score %>% 
  filter(X22q == -1) %>% 
  dplyr::select(X) %>%
  unlist(use.names = F) -> chr_22q_loss

df_aneu_score %>% 
  filter(X22q == 0) %>% 
  dplyr::select(X) %>% 
  unlist(use.names = F) -> chr_22q_disomy


# Primary setting
# based on the aneuploidy score
# 22q monosomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_22q_loss) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 25) %>%
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_22q_loss_pri_flt

# 22q disomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_22q_disomy) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 25) %>% 
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_22q_diso_pri_flt  


# RNAi dependency
RNAi_df_flt %>% 
  filter(X == "DDX5") -> DDX5_dep_RNAi

DDX5_dep_RNAi_t = as.data.frame(t(DDX5_dep_RNAi))

# rename the dep 
colnames(DDX5_dep_RNAi_t) = c("DDX5_dep")

DDX5_dep_RNAi_t$DepID = rownames(DDX5_dep_RNAi_t)
DDX5_dep_RNAi_t_edt = DDX5_dep_RNAi_t[-c(1),]


# assign the 22q loss and no loss 
DDX5_dep_RNAi_t_edt %>% 
  dplyr::mutate(DDX5_dep = as.numeric(DDX5_dep)) %>% 
  filter(DepID %in% c(chr_22q_loss_pri_flt, chr_22q_diso_pri_flt)) %>% 
  dplyr::mutate(condition = ifelse(DepID %in%chr_22q_loss_pri_flt ,"22q_loss", "22q_non_loss")) -> DDX5_dep

# plot the dependency 
DDX5_dep %>% 
  ggplot(aes(x = condition, y = DDX5_dep)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  theme(legend.position = "none") +
  labs(y = "DDX5_dep (RNAi)")

# p value =  2.26e-10
# only a few cell line in the RNAi
t.test(DDX5_dep[DDX5_dep$condition == "22q_loss",]$DDX5_dep, DDX5_dep[DDX5_dep$condition == "22q_non_loss",]$DDX5_dep)

# test for the correlation between 

# dependency and exp correlation analysis
DDX5_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  DDX5_dep %>% 
    rename(DepMap_ID = DepID) %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,1],comb_tmp[,4])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})


# combine the list
DDX5_exp_corr_cmb = do.call( rbind,DDX5_exp_corr)
colnames(DDX5_exp_corr_cmb) = c("gene", "p_value", "corr")


library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
DDX5_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(DDX5_exp_corr_cmb,gene%in% c("DDX17")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    #label.padding = 0.25,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
  ) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[2],"grey",colors[1])) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         legend.position = "none")



DDX17_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "DDX17"))]
DDX5_dep %>% 
  rename(DepMap_ID = DepID) %>% 
  left_join(DDX17_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(DDX5_dep)) %>% 
  filter(!is.na(DDX17)) -> cmb_df
cmb_df$pc = predict(prcomp(~DDX5_dep+DDX17, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = DDX17, y = DDX5_dep)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = 1)) +
  theme(legend.position = "none") +
  labs( y = "DDX5_dep (RNAi)", x = "DDX17 exp (TPM)") +
  geom_text( x = 5.8, y = -1.1, label = "R = 0.2804749  ")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value =3.69e-07", x = 5.8 ,y =-1.2 )
# correlation test
#p-value = 3.69e-07
cor_result = cor.test(cmb_df$DDX5_dep, cmb_df$DDX17)


######################### 5q loss ###########################
#get the 5q loss cells
df_aneu_score %>% 
  filter(X5q == -1) %>% 
  dplyr::select(X) %>%
  unlist(use.names = F) -> chr_5q_loss

df_aneu_score %>% 
  filter(X5q == 0) %>% 
  dplyr::select(X) %>% 
  unlist(use.names = F) -> chr_5q_disomy

# Primary setting
# based on the aneuploidy score
# 5q monosomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_5q_loss) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 26) %>%
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_5q_loss_pri_flt

# 5q disomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_5q_disomy) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 26) %>% 
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_5q_diso_pri_flt  


# RNAi dependency
RNAi_df_flt %>% 
  filter(X == "TYMS") -> TYMS_dep_RNAi

TYMS_dep_RNAi_t = as.data.frame(t(TYMS_dep_RNAi))

# rename the dep 
colnames(TYMS_dep_RNAi_t) = c("TYMS_dep")

TYMS_dep_RNAi_t$DepID = rownames(TYMS_dep_RNAi_t)
TYMS_dep_RNAi_t_edt = TYMS_dep_RNAi_t[-c(1),]


# assign the 22q loss and no loss 
TYMS_dep_RNAi_t_edt %>% 
  dplyr::mutate(TYMS_dep = as.numeric(TYMS_dep)) %>% 
  filter(DepID %in% c(chr_5q_loss_pri_flt, chr_5q_diso_pri_flt)) %>% 
  dplyr::mutate(condition = ifelse(DepID %in%chr_5q_loss_pri_flt ,"5q_loss", "5q_non_loss")) -> TYMS_dep

# plot the dependency 
TYMS_dep %>% 
  ggplot(aes(x = condition, y = TYMS_dep)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  theme(legend.position = "none") +
  labs(y = "TYMS_dep (RNAi)")

# p value =   0.1897
# only a few cell line in the RNAi
t.test(TYMS_dep[TYMS_dep$condition == "5q_loss",]$TYMS_dep, TYMS_dep[TYMS_dep$condition == "5q_non_loss",]$TYMS_dep)

# test for the correlation between 

# dependency and exp correlation analysis
TYMS_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  TYMS_dep %>% 
    rename(DepMap_ID = DepID) %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,1],comb_tmp[,4])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})


# combine the list
TYMS_exp_corr_cmb = do.call( rbind,TYMS_exp_corr)
colnames(TYMS_exp_corr_cmb) = c("gene", "p_value", "corr")


library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
TYMS_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(TYMS_exp_corr_cmb,gene%in% c("DHFR")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    #label.padding = 0.25,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
  ) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[2],"grey",colors[1])) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         legend.position = "none")



DHFR_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "DHFR"))]
TYMS_dep %>% 
  rename(DepMap_ID = DepID) %>% 
  left_join(DHFR_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(TYMS_dep)) %>% 
  filter(!is.na(DHFR)) -> cmb_df
cmb_df$pc = predict(prcomp(~TYMS_dep+DHFR, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = DHFR, y = TYMS_dep)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = 1)) +
  theme(legend.position = "none") +
  labs( y = "TYMS_dep (RNAi)", x = "DHFR exp (TPM)") +
  geom_text( x = 3, y = -2, label = "R = 0.2425942")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value =1.52e-05-07", x = 3 ,y =-2.2)
# correlation test
#p-value =  1.52e-05
cor_result = cor.test(cmb_df$TYMS_dep, cmb_df$DHFR)


# project: RNAi analysis for the paralog dependency on DDX3Y loss
# Author: Yi
# Date: 11.June.2025

# library
library(tidyverse)
library(MetBrewer)
library(ggthemes)
library(ggpubr)


## load the data

#RNAi dataset
RNAi_df = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /Depmap_RNAi/D2_combined_gene_dep_scores.csv")

# sample information
RNAi_sample_info = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /Depmap_RNAi/sample_info.csv")

# get the CCLE sample info
CCLE_sample_info = read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /Depmap_RNAi/CCLE_sample_info_file_2012-10-18.txt")

# get the Depmap sample info
Depmapsample_info = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /Depmap_RNAi/Model.csv")

#get the colname of RNAi
RNAi_colname = colnames(RNAi_df)

head(RNAi_colname)
# to get the cell line name and its corresponding Depmap id
cell_line_names = c()
for (i in c(2:length(RNAi_colname))) {
  # get the name 
  if(i <= 15){
    name = RNAi_colname[i]
    
    # split the name without "X"
    trimmed_name = str_sub(name, start = 2, end = nchar(name))
    
    cell_line_names = c(cell_line_names,trimmed_name )
  }else{
    trimmed_name = RNAi_colname[i]
    
    cell_line_names = c(cell_line_names,trimmed_name )
    
  }
  
}

# get the cell name to the depmap id
Depmapsample_info_ref = Depmapsample_info[,c("ModelID", "CCLEName")]

cell_line_name_df = as.data.frame(cell_line_names)
colnames(cell_line_name_df) = "CCLEName"
cell_line_name_df %>% 
  left_join(Depmapsample_info_ref, by = "CCLEName") -> join_tbl_CCLE_dep_name

# try to get rid of cell line without the info from Depmap
rm_col_index = as.integer(which(is.na(join_tbl_CCLE_dep_name$ModelID))) + 1

# remove the NA name
RNAi_df_flt = RNAi_df[,-rm_col_index]

# rename the colname into Depmap
colnames(RNAi_df_flt)[2:dim(RNAi_df_flt)[2]] = join_tbl_CCLE_dep_name$ModelID[!is.na(join_tbl_CCLE_dep_name$ModelID)]

# get the gene name and reform 
gene_name = RNAi_df_flt$X
gene_list = c()
for (i in c(1:length(gene_name))) {
  name = gene_name[i]
  
  gene = str_split(name, pattern = " ")[[1]][1]
  
  gene_list = c(gene_list, gene)
}

# rename the X column to the gene name
RNAi_df_flt$X = gene_list


# load the data from aneuploidy score

# dataset
# CCLE aneuploid status annotation
aneuploid_dep = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/aneuploidy_scores.csv")
# CCLE cell line with each chromsome arm gain or loss info
df_aneu_score = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/arm_call_scores.csv")
#cell line annotation
cell_line_id = read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/Cell_lines_annotations_20181226.txt")

# color set
colors = c("#EF8536","#3A76AF")


# correlation dependency
## analysis here: correlation and dep and exp
# retrieve gene expression data
CCLE_exp = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")

# change the name of CCL2_exp data
colnames(CCLE_exp)[1] = "DepMap_ID"
#dim(df_gene_dep)

# gene the gene name 
#dim(df_gene_dep)
genelist_exp = sub("\\..*","",colnames(CCLE_exp)[2:dim(CCLE_exp)[2]])

# replace gene to colnames
colnames(CCLE_exp)[2:dim(CCLE_exp)[2]] = genelist_exp


# Look at DDX3X loss only and the dataset is missing, we need to redefine the loss of chromosome Y loss

# gene copy number
CCLE_gene_cnv = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /PortalOmicsCNGeneLog2.csv")

# cell and patient info
cell_line_id = read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/Cell_lines_annotations_20181226.txt")

# rename the gene name on CCLE_gene_CNV
# change the name of data
colnames(CCLE_gene_cnv)[1] = "DepMap_ID"
#dim(df_gene_dep)

# gene the gene name 
dim(CCLE_gene_cnv)
genelist_Dep = sub("\\..*","",colnames(CCLE_gene_cnv)[2:dim(CCLE_gene_cnv)[2]])

# replace gene to colnames
colnames(CCLE_gene_cnv)[2:dim(CCLE_gene_cnv)[2]] = genelist_Dep
head(CCLE_gene_cnv[,1:6])

library(biomaRt)
# retrieve Y chromosome genes
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
results <- getBM(attributes = c("chromosome_name", "entrezgene_id", "hgnc_symbol"),
                 filters = "chromosome_name", values = "Y", mart = mart)

Y_chr_genename = results$hgnc_symbol[results$hgnc_symbol != ""]

# get the table as Y chromosome cnv
Y_cnv = CCLE_gene_cnv[,c(1,which(colnames(CCLE_gene_cnv) %in% Y_chr_genename))]

# get the male cell line
cell_line_id %>% 
  dplyr::filter(Gender != "female") %>% 
  dplyr::select(depMapID) %>% 
  unlist(use.names = F) -> male_cellline_flt

# get the male Y chromomr cnv df
Y_cnv %>% 
  dplyr::filter(DepMap_ID %in% male_cellline_flt) -> Y_cnv_male

# Y chr loss sample
Y_loss_sample = Y_cnv_male[rowSums(is.na(Y_cnv_male))>25,]$DepMap_ID 
# Y chr normal sample
Y_normal = Y_cnv_male$DepMap_ID[!(Y_cnv_male$DepMap_ID %in% Y_loss_sample)]


# The gene we found in Y chr loss dependency is DDX3X and 
RNAi_df_flt %>% 
  filter(X == "DDX3X") -> DDX3X_dep_RNAi

DDX3X_dep_RNAi_t = as.data.frame(t(DDX3X_dep_RNAi))

# rename the dep 
colnames(DDX3X_dep_RNAi_t) = c("DDX3X_dep")

DDX3X_dep_RNAi_t$DepID = rownames(DDX3X_dep_RNAi_t)
DDX3X_dep_RNAi_t_edt = DDX3X_dep_RNAi_t[-c(1),]

# assign the 8p loss and no loss 
DDX3X_dep_RNAi_t_edt %>% 
  dplyr::mutate(DDX3X_dep = as.numeric(DDX3X_dep)) %>% 
  filter(DepID %in% c(Y_loss_sample, Y_normal)) %>% 
  dplyr::mutate(condition = ifelse(DepID %in%Y_loss_sample ,"Y_loss", "Y_non_loss")) -> DDX3X_dep

# plot the dependency 
DDX3X_dep %>% 
  ggplot(aes(x = condition, y = DDX3X_dep)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  theme(legend.position = "none") +
  labs(y = "DDX3X_dep (RNAi)")


# only a few cell line in the RNAi
# p = 0.04794
t.test(DDX3X_dep[DDX3X_dep$condition == "Y_loss",]$DDX3X_dep, DDX3X_dep[DDX3X_dep$condition == "Y_non_loss",]$DDX3X_dep)

# test for the correlation between 

# dependency and exp correlation analysis
DDX3X_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  DDX3X_dep %>% 
    rename(DepMap_ID = DepID) %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,1],comb_tmp[,4])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})


# combine the list
DDX3X_exp_corr_cmb = do.call( rbind,DDX3X_exp_corr)
colnames(DDX3X_exp_corr_cmb) = c("gene", "p_value", "corr")


library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
DDX3X_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(DDX3X_exp_corr_cmb,gene%in% c("DDX3Y")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    #label.padding = 0.25,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
  ) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[2],"grey",colors[1])) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         legend.position = "none")


# correlation plot
# based on the expression and correlation of DDX3Y

DDX3Y_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "DDX3Y"))]
DDX3X_dep %>% 
  rename(DepMap_ID = DepID) %>% 
  left_join(DDX3Y_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(DDX3X_dep)) %>% 
  filter(!is.na(DDX3Y)) -> cmb_df
cmb_df$pc = predict(prcomp(~DDX3X_dep+DDX3Y, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = DDX3Y, y = DDX3X_dep)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = -1)) +
  theme(legend.position = "none") +
  labs( y = "DDX3X_dep (RNAi)", x = "DDX3Y exp (TPM)") +
  geom_text( x = 2, y = 0.58 , label = "R = 0.3286202 ")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value = 2.01e-08", x = 2, y =0.5) 

# correlation test
#p-value = 2.01e-08
cor_result = cor.test(cmb_df$DDX3X_dep, cmb_df$DDX3Y)


