
#--------------------------------------------------------------------------------------------------------------------
#
# Extract proteins using the search terms for the proteins involved in endothelial activation and inflammation 
#
# By Zihan Sun 
#
#--------------------------------------------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
#rm(list=ls())

setwd("~/immune_search_terms")

# grep returns the index of the pattern 
# grepl returns the true or false of whether the pattern exists in the list 

# load the keyword list (organized from Table S1)
key_list = fread("Keywords_list_04032024.txt")

key_list[key_list == ""] = NA 

key_list$Inflammation = str_replace(key_list$Inflammation, ' ""',"")
key_list$Inflammation = str_replace(key_list$Inflammation, '""',"")
Inflammation = unique(na.omit(key_list$Inflammation))
#print.noquote(Inflammation)
cat(paste(Inflammation, ","))

key_list$BBB = gsub(' ""',"", key_list$BBB)
key_list$BBB = gsub('""',"", key_list$BBB)
BBB = unique(na.omit(key_list$BBB))
cat(paste(BBB, ","))

key_list$`Endothelial dysfunction` = gsub(' ""',"", key_list$`Endothelial dysfunction`)
key_list$`Endothelial dysfunction` = gsub('""',"", key_list$`Endothelial dysfunction`)
Endo = unique(na.omit(key_list$`Endothelial dysfunction`))
cat(paste(Endo, ","))

key_list$`Oxidative stress` = gsub(' ""',"", key_list$`Oxidative stress`)
key_list$`Oxidative stress` = gsub('""',"", key_list$`Oxidative stress`)
Oxi = unique(na.omit(key_list$`Oxidative stress`))
cat(paste(Oxi, ","))

key_list$`Neuronal and glial function` = gsub(' ""',"", key_list$`Neuronal and glial function`)
key_list$`Neuronal and glial function` = gsub('""',"", key_list$`Neuronal and glial function`)
Neuro = unique(na.omit(key_list$`Neuronal and glial function`))
cat(paste(Neuro , ","))

key_list$`Clotting, vascular remodeling, and metabolic reprogramming ` = gsub(' ""',"", key_list$`Clotting, vascular remodeling, and metabolic reprogramming `)
key_list$`Clotting, vascular remodeling, and metabolic reprogramming ` = gsub('""',"", key_list$`Clotting, vascular remodeling, and metabolic reprogramming `)
Vas = unique(na.omit(key_list$`Clotting, vascular remodeling, and metabolic reprogramming `))
cat(paste(Vas, ","))

exposure = c(Inflammation, BBB, Endo, Oxi, Neuro, Vas)

#----------------------------------------------
# 
# Get the related protein list 
#
#----------------------------------------------

library(readxl)

#=========
# SOMA
#=========

# this sheet is the Supplemental Table 2 from the Eldjarn et al. Nature 2023 paper
soma_info = read_excel("~/OneDrive - University of Cambridge/Summary Sheets/CSF_brain_plasma_protein_list/plasma_pQTL_sumstats/Olink_Soma_Comparison/Olink_Soma_Comparison.xlsx", sheet = "ST2_somascan_assays")

colnames(soma_info) = soma_info[2, ]
soma_info = soma_info[c(-1, -2), ]

immune_related_soma  <- subset(soma_info, grepl(regex(paste(exposure, collapse = "|"), ignore_case = T), target_full_name, ignore.case=T)
                               | grepl(regex(paste(exposure, collapse = "|"), ignore_case = T), target_name, ignore.case=T)
                               | grepl(regex(paste(exposure, collapse = "|"), ignore_case = T), gene_name, ignore.case=T))

immune_related_soma  <- immune_related_soma %>% 
  select("seqid", "uniprot", "gene_name", "target_name", "target_full_name", "norm_cv") %>%
  rename(target_name_soma = target_name, cv_soma = norm_cv) %>%
  mutate(uniprot = str_replace(uniprot, ", ", "")) %>%
  mutate(uniprot = str_replace(uniprot, ",", "")) %>%
  mutate(uniprot = str_replace(uniprot, ",", "")) %>%
  mutate(uniprot = str_replace(uniprot, " ", "")) %>%
  mutate(uniprot = str_replace(uniprot, " ", "")) %>%
  mutate(uniprot = str_replace(uniprot, "  ", "")) 
# re-organize the assays which have multiple uniprot IDs 

immune_related_soma_multi  <- immune_related_soma %>% 
  filter(nchar(uniprot) > 6)


length(unique(immune_related_soma$seqid))
# 2144
# 2144 unique ASSAYS identified from keyword search from SomaScan Icelander 36K study
length(unique(immune_related_soma$uniprot))
# 2013 !!!

#=========
# OLINK
#=========

# this sheet is the Olink info sheet from the Supplemental Materials of the Eldjarn et al. Nature 2023 paper
olink_info = fread("~/OneDrive - University of Cambridge/Summary Sheets/CSF_brain_plasma_protein_list/plasma_pQTL_sumstats/Olink_Soma_Comparison/olink_protein_info.txt")

olink_info = olink_info %>%
  mutate(uniprot = str_replace(uniprot, "_", "")) %>%
  mutate(uniprot = str_replace(uniprot, "_", ""))

olink_info_multi  <- olink_info %>% 
  filter(nchar(uniprot) > 6)

#protein_full_names = olink_info[, 1:4]

immune_related_olink  <- subset(olink_info, grepl(regex(paste(exposure, collapse = "|"), ignore_case = T), target_full_name, ignore.case=T)
                                | grepl(regex(paste(exposure, collapse = "|"), ignore_case = T), gene_name, ignore.case=T))

immune_related_olink  <- immune_related_olink %>% 
  select("oid", "uniprot", "gene_name", "target_full_name", "cv") %>%
  rename(cv_olink = cv) 

immune_related_olink_multi = immune_related_olink %>% 
  filter(nchar(uniprot) > 6)

length(unique(immune_related_olink$oid))
# 1364
# 1364 unique assays identified from keyword search from UKB-PPP all panels
length(unique(immune_related_olink$uniprot))
# 1355 !!!

#-----------------------------------
# 
# plus Olink Inflam Panel I and II 
#
#-----------------------------------

ukb_ppp = fread("~/OneDrive - University of Cambridge/Summary Sheets/CSF_brain_plasma_protein_list/plasma_pQTL_sumstats/UKB-PPP/plasma_UKB-PPP_pQTL_list.txt")
length(unique(ukb_ppp$UniProt))
# 2922 unique proteins assessed in UKB-PPP which had pQTL

panel_info = ukb_ppp %>%
  select(`Olink ID`, `Protein panel`, UniProt) %>%
  mutate(UniProt = str_replace(UniProt, ";", "")) %>%
  mutate(UniProt = str_replace(UniProt, ";", ""))

# 2941 analytes

colnames(panel_info) = c("oid", "panel", "uniprot")

olink_info_panel = left_join(panel_info, olink_info, by = "oid")

olink_inflam = olink_info_panel %>% 
  filter(panel == "Inflammation" | panel == "Inflammation_II")
length(unique(olink_inflam$uniprot.y))
# 736 unique assays in Inflammation Panels I and II

immune_related_olink_inflam  <- olink_inflam %>% 
  select("oid", "uniprot.y", "gene_name", "target_full_name", "cv") %>%
  rename(cv_olink = cv, uniprot = uniprot.y) 

immune_related_olink_inflam_multi = immune_related_olink_inflam %>% 
  filter(nchar(uniprot) > 6)

length(unique(immune_related_olink_inflam$uniprot))
# 736 !!!!
# 736 unique assays in Olink Inflam Panel I and II

length(unique(immune_related_olink_inflam$oid))

#===========================================
#
# Merge Soma + Olink + Olink_Inflam
#
#===========================================

immune_related = full_join(immune_related_olink, immune_related_olink_inflam, by = c("uniprot", "gene_name", "target_full_name", "cv_olink", "oid"))
immune_related = full_join(immune_related, immune_related_soma, by = c("uniprot"))

immune_related = immune_related %>%
  select(uniprot, oid, seqid, gene_name.x, target_full_name.x, cv_olink, cv_soma) %>%
  mutate(cv_soma = as.numeric(cv_soma))

length(unique(immune_related$uniprot))
# 2812
# 2812 unique proteins identified as inflammation or BBB related 

#========================================================================================
# (check how many duplicates are removed)

overlap_s1 = inner_join(immune_related_olink, immune_related_olink_inflam, by = c("uniprot", "gene_name", "target_full_name", "cv_olink", "oid"))
length(unique(overlap_s1$uniprot))
# 393 unique proteins are overlapped in olink-Inflam panel and olink key word search

immune_related_olink_both_sources = full_join(immune_related_olink, immune_related_olink_inflam, by = c("uniprot", "gene_name", "target_full_name", "cv_olink", "oid"))
overlap_s2 = inner_join(immune_related_olink_both_sources, immune_related_soma, by = c("uniprot"))
length(unique(overlap_s2$uniprot))
# 899 unique proteins are overlapped in olink panel and soma key word search

#overlap_total = full_join(overlap_s1, overlap_s2, by = c("uniprot", "gene_name", "target_full_name", "cv_olink", "oid") )
#length(unique(overlap_total$uniprot))
# 948 unique proteins are overlapped either two or all three of the search methods

#========================================================================================
# for the overlapped proteins between olink and soma, use the results from olink
# get the list of proteins available in Olink

immune_related_keep_olink = immune_related %>% 
  filter(!is.na(oid)) %>%
  group_by(uniprot) %>% 
  arrange(cv_olink) %>% 
  ungroup() %>%
  distinct(uniprot, .keep_all = T)

length(unique(immune_related_keep_olink$uniprot))
# 1698

# get the list of proteins unique to Soma v4

immune_related_uniq_soma = immune_related %>% 
  filter(is.na(oid)) %>%
  group_by(uniprot) %>% 
  arrange(cv_soma) %>% 
  ungroup() %>%
  distinct(uniprot, .keep_all = T)

length(unique(immune_related_uniq_soma$uniprot))
# 1114

###############################################
# 
# Perform cis-pQTL, CV, and LOD filtering  
#
###############################################

#------------------------------------------
# UKB - prepare UKB pQTL and NPX QC info
# check ST info and retain relevant columns 
#------------------------------------------------

# 1. primary sig cis signals

ukb_ppp_replicate_cis = fread("~/OneDrive - University of Cambridge/Summary Sheets/CSF_brain_plasma_protein_list/plasma_pQTL_sumstats/UKB-PPP/UKB-PPP_replicate_pQTL.txt")
# 1955 primary cis associations

length(unique(ukb_ppp_replicate_cis$`Assay Target`))
# 1954

ukb_ppp_replicate_cis[which(duplicated((ukb_ppp_replicate_cis$`Assay Target`))), ]
#EBI3_IL27 have two UniProt IDs

length(unique(ukb_ppp_replicate_cis$`Target UniProt`))
# 1953

ukb_ppp_replicate_cis[which(duplicated((ukb_ppp_replicate_cis$`Target UniProt`))), ]
#IL12A_IL12B, IL12B have the same uniprot id
#NPPB, NTproBNP have the same uniprot id

#------------------------------------------------------------
# 2. primary sig cis and trans signals and external validation

ukb_ppp_discover_external_validation_signals = fread("~/OneDrive - University of Cambridge/Summary Sheets/CSF_brain_plasma_protein_list/plasma_pQTL_sumstats/UKB-PPP/UKB-PPP_ST15_discover_external_validation_signals.txt")
# 14287 in total (matched with the result reported in the Sun et al. paper)

ukb_ppp_discover_external_validation_signals = ukb_ppp_discover_external_validation_signals %>% 
  separate(`UKBPPP ProteinID`, c("gene_name", "uniprot", "oid", "v1"), sep = ":") %>%
  select("CHROM", "GENPOS (hg38)",  "oid", "uniprot", "gene_name", "MHC", "novel/replicated", "cis/trans") %>%
  rename(external_validation = `novel/replicated`)

#------------------------------------------------------------
# 3. primary sig signals and internal validation

ukb_ppp_discover_internal_validation_signals = fread("~/OneDrive - University of Cambridge/Summary Sheets/CSF_brain_plasma_protein_list/plasma_pQTL_sumstats/UKB-PPP/UKB-PPP_ST9_discover_replicate_signals.txt")

# 1.3 * 10^(−5) is the p-threshold used for replication; 1.7 * 10^(-11) is used for discovery
neg_log10p_dis = -log10(1.7 * 10^(-11))
neg_log10p_rep = -log10(1.3 * 10^(-5))
#neg_log10p = 4.886057
#neg_log10p = 4.89

#10^(-4.886057) = 1.299999e-05
#10^(-10.8) = 1.299999e-05

ukb_ppp_discover_internal_validation_signals = ukb_ppp_discover_internal_validation_signals %>% 
  filter(!is.na(`BETA (discovery, wrt. A1)`)) %>%
  separate(`UKBPPP ProteinID`, c("gene_name", "uniprot", "oid", "v1"), sep = ":") %>%
  filter(`BETA (discovery, wrt. A1)` != 0 & `BETA (replication)` != 0) %>%  # no beta equals 0
  mutate(beta_dir = case_when(`BETA (discovery, wrt. A1)` < 0 & `BETA (replication)` < 0 ~ "neg",
                              `BETA (discovery, wrt. A1)` > 0 & `BETA (replication)` > 0 ~ "pos",
                              `BETA (discovery, wrt. A1)` < 0 & `BETA (replication)` > 0 ~ "inconsistent",
                              `BETA (discovery, wrt. A1)` > 0 & `BETA (replication)` < 0 ~ "inconsistent")) %>%
  mutate(`novel/replicated` = case_when (beta_dir == "neg" & `log10(p) (replication)` > neg_log10p_rep ~ "replicated",  #!!
                                         beta_dir == "pos" & `log10(p) (replication)` >  neg_log10p_rep ~ "replicated",
                                         beta_dir == "neg" & `log10(p) (replication)` <=  neg_log10p_rep ~ "novel",  #!!
                                         beta_dir == "pos" & `log10(p) (replication)` <=  neg_log10p_rep ~ "novel",
                                         beta_dir == "inconsistent" ~ "novel")) %>%
  select("CHROM", "GENPOS (hg38)",  "oid", "uniprot", "gene_name", "MHC", "novel/replicated", "cis/trans", "log10(p) (discovery)") %>%
  rename(internal_validation = `novel/replicated`)

table(ukb_ppp_discover_internal_validation_signals$internal_validation)

#novel replicated 
#4455       9832 

4455 + 9832 = 14287
1955 + 12332 = 14287
  
1869 + 7906 = 9775

# note that the numbers that we got is slightly higher than the numbers reported in the paper 9832 vs 9775 Num. of associations are replicated 
# This is possibly because of the rounding of the -log10(p-value)
# however, we were unable to recover the exact p-values based on the summary data that the authors reported 
# since our analysis is meant to be inclusive at this step, we used the replicated list found in our curation based on summary info in ST9. 

#----------------------------------------------------------------------------
# get the full replication list (internal + external validation)
#----------------------------------------------------------------------------

ukb_ppp_signals = inner_join(ukb_ppp_discover_internal_validation_signals,
                             ukb_ppp_discover_external_validation_signals, 
                             by = c("CHROM", "GENPOS (hg38)",  "oid", "uniprot", "gene_name", "MHC", "cis/trans") )

# all 14,287 signals are matched up! 

# as long as the signal is replicated either in the replication set in ukb or in an external dataset, they are counted as replicated
ukb_ppp_signals = ukb_ppp_signals %>% 
  mutate(replication = ifelse(internal_validation == "replicated" | external_validation == "replicated", 
                              "yes", "no"))

#---------------------
# check Olink LOD 
# % data fail QC
#---------------------

library(readxl)
olink_lod = read_excel("~/OneDrive - University of Cambridge/Summary Sheets/CSF_brain_plasma_protein_list/plasma_pQTL_sumstats/Olink_Soma_Comparison/Olink_Soma_Comparison.xlsx", sheet = "ST24_olink_perc_below_lod")
olink_lod = olink_lod[c(-1, -2), ]
colnames(olink_lod) = c("oid", "percent_below_LOD")
olink_lod$percent_below_LOD = as.numeric(olink_lod$percent_below_LOD)

#--------------
# concise info
#--------------

ukb_ppp_sig_var_info = ukb_ppp_signals %>% 
  inner_join(olink_lod, by = "oid") %>%
  inner_join(olink_info, by = c("oid")) %>%
  rename(uniprot = uniprot.y, gene_name = gene_name.y) %>% 
  select(-c("uniprot.x", "gene_name.x"))   # keep the gene_name and uniprot id from olink_info file


length(unique(ukb_ppp_sig_var_info$uniprot))
# primary singals (cis / trans) were identified for 2414 unique proteins 


#---------------------------------------------------------
# perform step-wise selection in ukb-ppp proteins
#---------------------------------------------------------

ukb_immune_list = immune_related_keep_olink %>%      
  select("oid") %>% 
  inner_join(ukb_ppp_sig_var_info, by = "oid")

length(unique(immune_related_keep_olink$oid)) # 1698 unique olink assays are immune-related
length(unique(ukb_immune_list$oid))  #primary associations were identified in the UKB-discovery cohort for 1448 assays

length(unique(immune_related_keep_olink$uniprot)) # 1698 unique olink assays are immune-related
length(unique(ukb_immune_list$uniprot))       # 1448 unique proteins 


ukb_immune_list_cis = ukb_immune_list %>%      
  filter(`cis/trans` == "cis")
length(unique(ukb_immune_list_cis$uniprot))  # 1222 proteins have cis-associations identified 

ukb_immune_list_cis_rep = ukb_immune_list_cis %>%      
  filter(replication == "yes")
length(unique(ukb_immune_list_cis_rep$uniprot))  # 1183 proteins have been replicated either internally within the UKB-pop or externally in other cohort

ukb_immune_list_cis_rep_nomhc = ukb_immune_list_cis_rep %>%      
  filter(MHC==0)
length(unique(ukb_immune_list_cis_rep_nomhc$uniprot))  # 1162 proteins whose replicated cis-pQTL does not occur in the MHC region

ukb_immune_list_cis_rep_nomhc_cv = ukb_immune_list_cis_rep_nomhc %>%      
  filter(cv < 20)
length(unique(ukb_immune_list_cis_rep_nomhc_cv$uniprot)) # 1119 proteins have cv <20

ukb_immune_list_cis_rep_nomhc_cv_lod = ukb_immune_list_cis_rep_nomhc_cv %>%      
  filter(percent_below_LOD < 0.05)
length(unique(ukb_immune_list_cis_rep_nomhc_cv_lod$uniprot)) # 979 proteins have >= 5% above LOD 

ukb_immune_list_cis_rep_nomhc_cv_lod_Xchr = ukb_immune_list_cis_rep_nomhc_cv_lod %>%      
  filter(CHROM != 23)
length(unique(ukb_immune_list_cis_rep_nomhc_cv_lod_Xchr$uniprot)) # 956 proteins have >= 5% above LOD --> 23 proteins were excluded 

# 7 proteins CD63, FLT1, IL12A_IL12B, NBN, PPBP, PTRHD1, UBXN1 were additionally excluded during the identification of MR instruments.
# these 7 proteins do not have any instrument available when harmonized with any of the outcome datasets 
# the code for finding the MR instruments is documented the MR analysis code

ukb_immune_list_cis_rep_nomhc_cv_lod_Xchr = ukb_immune_list_cis_rep_nomhc_cv_lod_Xchr %>%      
  filter(gene_name != "CD63" & gene_name != "FLT1" & gene_name != "IL12A_IL12B" & 
           gene_name != "NBN" & gene_name != "PPBP" & gene_name != "PTRHD1" & gene_name != "UBXN1")
length(unique(ukb_immune_list_cis_rep_nomhc_cv_lod_Xchr$uniprot)) # 949

#================================
# DONE WITH UKB DATA PROCESSING 
# Cheers !!!!!
#================================

#----------------------------------------------------
#
# Filter the candidate proteins unique to SomaScan v4 
#
#----------------------------------------------------

#----------
# check 
#----------

icelander = fread("~/OneDrive - University of Cambridge/Summary Sheets/CSF_brain_plasma_protein_list/plasma_pQTL_sumstats/Icelander/Icelander36K_pQTL_SMPnorm.txt")
length(unique(icelander$SeqId))
# 4809 assays were included 

length(unique(icelander$UniProt))
# 4576 with a unique UniProt n ID has either >0 cis or trans pQTLs (or both) identified

icelander_sentinel = icelander %>% filter(`rank
locus` == 1)
# 24,736 sentinel pQTL in total 

icelander_sentinel_cis = icelander_sentinel %>% filter(`cis
trans`== "cis") 
# 2120 cis sentinel pQTLs --> number matched with their paper

icelander_sentinel_trans = icelander_sentinel %>% filter(`cis
trans`== "trans") 
# 22,616 trans sentinel pQTLs --> number matched with their paper

colnames(icelander) = c("locus_ID_global", "locus_ID_prot", "gene_prot", "chr_prot", "strand_prot", "TSS_prot", "shortname_prot", "long_name_prot", "UniProt", 
                        "variant", "LD_class_size", "LD_class", "chr_var", "pos_var", "region_start", "region_end", 
                        "Amin", "Amaj", "info", "MAF_PC", "cis_trans", "distTSS", "closest_genes", "cis_eqtl_coding_genes", "beta_adj", 
                        "beta_unadj", "mLog10pval_adj", "mLog10pval_unadj", "rank_locus", "count_locus", "count_cis_eqtl_same_gene", 
                        "prop_same_direction_cis_eqtl_same_gene", "count_cis_eqtl_diff_gene", "any_coding_same_gene", "any_coding_diff_gene", 
                        "rank_proteins_per_pqtl", "count_proteins_per_pqtl", "count_cis_proteins_per_pqtl", "rank_pqtl_per_protein", 
                        "count_pqtl_per_protein", "reported", "reported_rsid", "reported_p", "reported_pmid", "reported_r2", "SeqId", 
                        "r2_joint_pqtl", "effect_in_joint_model", "p-value_joint_model")


#------------------------------
# Find MHC yes / no
#------------------------------

# chromosome 6: 25.5–34.0 Mb

mhc_soma_marker = icelander %>% 
  filter(chr_var == "chr6") %>% 
  filter(pos_var >= 25500000 & pos_var <= 34000000) %>%
  filter(cis_trans == "cis") %>%
  filter(reported == "Y") 

mhc_id = unique(mhc_soma_marker$SeqId)
# 8 assays are in this list

#---------------------------------------------------------
# perform step-wise selection in Soma proteins
#---------------------------------------------------------

length(unique(immune_related_uniq_soma$seqid))  #1114 unique soma assays related to immunity

length(unique(immune_related_uniq_soma$uniprot))  #1114 unique soma proteins related to immunity

immune_related_uniq_soma= immune_related_uniq_soma %>% 
  dplyr::mutate(SeqId = stringr::str_replace(seqid, "SeqId.", "")) %>%
  dplyr::mutate(SeqId = stringr::str_replace(SeqId, "-", "_")) %>%
  select(SeqId, cv_soma)

soma_immune_list = inner_join(immune_related_uniq_soma, icelander, by = "SeqId")

length(unique(soma_immune_list$SeqId))  #primary associations were identified in the icelander cohort for 1051 assays
length(unique(soma_immune_list$UniProt))       # 1051 Uniprot IDs

soma_immune_list_cis = soma_immune_list %>%      
  filter(cis_trans == "cis")
length(unique(soma_immune_list_cis$UniProt))  # 356 proteins have cis-associations identified 

soma_immune_list_cis_rep = soma_immune_list_cis %>%      
  filter(reported == "Y")
length(unique(soma_immune_list_cis_rep$UniProt))  # 53 proteins have been replicated externally in other cohort

soma_immune_list_cis_rep_nomhc = soma_immune_list_cis_rep %>%      
  filter(! SeqId %in% mhc_id)
length(unique(soma_immune_list_cis_rep_nomhc$UniProt))  # 49 proteins whose replicated cis-pQTL does not occur in the MHC region

soma_immune_list_cis_rep_nomhc_cv = soma_immune_list_cis_rep_nomhc %>%      
  filter(cv_soma < 20)
length(unique(soma_immune_list_cis_rep_nomhc_cv$UniProt)) # 47 proteins have cv <20


# !!!! soma do not have LOD issues

#================================
# get the final list of proteins
#================================

#--------
# Soma 
#-------

####################################################################################################################

soma_final= soma_immune_list_cis_rep_nomhc_cv %>% 
  select(SeqId, gene_prot, chr_prot, TSS_prot, shortname_prot, long_name_prot, UniProt) %>% 
  distinct(SeqId, .keep_all = T)

write.table(soma_final, "soma_final_list_47_05032024.txt", col.names=T, row.names=F, quote=F, sep="\t", append = F)


#--------
# Olink 
#--------

ukb_ppp = ukb_ppp %>% 
  rename(oid = `Olink ID`)

# get the format matched with the GWAS sumstats naming conventions

ukb_final= ukb_immune_list_cis_rep_nomhc_cv_lod_Xchr %>%
  select(gene_name, uniprot, oid) %>% 
  distinct(oid, .keep_all = T) %>%
  #mutate(`UKBPPP ProteinID` = paste(gene_name, uniprot, oid, "v1", sep = ":")) %>%
  left_join(ukb_ppp, by = "oid") %>%
  select(`UKBPPP ProteinID`, `Protein panel`) %>%
  mutate(sumstats_name = str_replace_all(`UKBPPP ProteinID`, ":", "_")) %>%
  mutate(sumstats_name = paste(sumstats_name, `Protein panel`, sep ="_")) 
#%>%select(sumstats_name)

write.table(ukb_final, "ukb_list_949.txt", col.names=T, row.names=F, quote=F, sep="\t", append = F)

#ukb_final[which(is.na(ukb_final)), ]
# all matched up (no NAs)


#==================================
# get the tss for each protein !!!!
#==================================

#ukb_tss = fread("ukb_gene_start_end.txt")
# this file is no longer needed

ukb_tss= ukb_immune_list_cis_rep_nomhc_cv_lod_Xchr %>%
  select(gene_name, uniprot, oid) %>% 
  distinct(oid, .keep_all = T) %>%
  left_join(ukb_ppp, by = "oid") %>%
  select(`UKBPPP ProteinID`, `Protein panel`, oid, `Gene CHROM`, `Gene start`, `Gene end`) %>%
  mutate(folder_name = str_replace_all(`UKBPPP ProteinID`, ":", "_")) %>%
  mutate(folder_name = paste(folder_name, `Protein panel`, sep ="_")) %>%
  mutate(file_name = paste(`UKBPPP ProteinID`, `Protein panel`, sep =":")) %>%
  select(folder_name, file_name, oid, `Gene CHROM`, `Gene start`, `Gene end`)

colnames(ukb_tss) = c("folder_name", "file_name", "oid", "CHR", "gene_start", "gene_end")

write.table(ukb_tss, "ukb_prot_tss_949.txt", col.names=T, row.names=F, quote=F, sep="\t", append = F)

table(ukb_tss$CHR)
# 4 assays include >1 uniprot proteins

