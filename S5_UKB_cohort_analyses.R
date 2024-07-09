
#-------------------------------------------------------------
# 
# Survival and cross-sectional analyses in the UKB-PPP cohort
# 
#-------------------------------------------------------------

library(data.table)
library(dplyr)
library(tidyr)

#rm(list=ls())

ukb_rdata = readRDS("datasets/UK_Biobank/phenotype/UK_Biobank_Olink/olink_data_reformatted_081123.rds")

dim(ukb_rdata)
# 55323  2925

length(unique(ukb_rdata$eid))  # 53070

table(ukb_rdata$ins_index)
# 0     2     3 
# 53026  1173  1124 
# 53026+1173+1124 = 55323

# get the baseline measure 
ukb_bl = ukb_rdata[which(ukb_rdata$ins_index == 0), ]  # 53026
dim(ukb_bl)
#53026  2925

bl_eid = as.data.frame(ukb_bl[ , "eid"])

fwrite(bl_eid, 
       "UKB_NOV2023/Olink_BL_participants_50K_eid.txt.gz", 
       col.names=T, row.names=F, quote=F, sep="\t", append = F)

#--------------------------
# continue qc baseline data
#--------------------------

# Step 1: 
# check the number of missing proteins per person
# margin =1 means across rows 

ukb_bl$N_missing_prot = apply(X = is.na(ukb_bl), MARGIN = 1, FUN = sum)
# get % missingness
ukb_bl$P_missing_prot = (ukb_bl$N_missing_prot / 2923 )*100

summary(ukb_bl$N_missing_prot)
summary(ukb_bl$P_missing_prot)

# histogram of the % missingness of protein measures across individuals at baseline
p1= hist(ukb_bl$P_missing_prot, main = "", xlab = "The percent of missing proteins per person",
         ylab = "No. of Individuals")

png(filename="UKB_NOV2023/hist_P_missing_prot.png")
plot(p1)
dev.off()

# Step 2: 
# check the number of missing persons per protein
# margin = 2 means across columns
N_missing_person = apply(X = is.na(ukb_bl[,2:1139]), MARGIN = 2, FUN = sum)
P_missing_person = (N_missing_person / 53026 )*100

#which(is.na(P_missing_person))
summary(N_missing_person)
summary(P_missing_person)

# histogram of the % sample missingness by proteins measured at baseline
p2 = hist(P_missing_person, main = "", xlab = " % sample missingness by proteins measured at baseline", 
          ylab = "No. of Proteins")

png(filename="UKB_NOV2023/hist_P_missing_person.png")
plot(p2)
dev.off()

# For the current analysis (21-02-2024), all the participants (n=53026) and all proteins (N=2923) are included in the analysis !!
# After merging with the phenotype datasets (downloaded and processed in Nov-2023), n = 53022 participants (ppts) remain. 

#===================================================================================================================

dat = readRDS("datasets/UK_Biobank/phenotype/UK_Biobank_NOV2023/UKB_nov2023_676165_clean.rds")
#colnames(dat)

# extract the relevant data fields from the full data sheet
df_use = dat[ , c("f.eid", "date_extraction", "date_baseline", 
                  
                  "age_baseline", "dob_year", "dob_month", "sex", "edu_yrs", "edu_cat_detailed", "edu_cat", "townsend",
                  
                  "bmi", "sbp", "sbp_sd", "dbp", "htn", "smallstat", "smallbin", "alcallstat", "alcallbin", 
                  
                  "bpmed", "cholmed", "insulinmed",
                  
                  "tchol", "tchol_mgdl", "hdl", "hdl_mgdl", "tchol_hdl_ratio", "ldl", "ldl_mgdl", 
                  "trig", "lntrig", "trig_mgdl", "lntrig_mgdl", "crp", "lncrp", "crp_mgdl", "lncrp_mgdl", "cog_global", "diabetes_t2d", "diabetes_t2d_selfreport",
                  
                  "fram_points", "fram_cvd_cat", "fram_cvd_cat3", 
                  
                  "rs429358_allele", "rs7412_allele", "apoe_isoform_detailed", "apoe_isoform", "apoe_e4_carrier",   
                  
                  "duration1", 
                  
                  "cog_reaction_time",                                
                  "cog_numeric_memory",                              
                  "cog_fluid_intelligence",                         
                  "cog_pairs_matching",                              
                  
                  "alg_allcause_dementia_date", 
                  "alg_alzheimer_date", 
                  "alg_vascular_dementia_date", 
                  "alg_stroke_date", 
                  "alg_ischstroke_date", 
                  
                  "fup_allcause_dementia_yrs_alg",                 
                  #"fup_allcause_dementia_yrs_man",
                  "ep1_allcause_dementia_alg",                      
                  #"ep1_allcause_dementia_man", 
                  
                  "fup_alzheimer_yrs_alg",                           
                  #"fup_alzheimer_yrs_man", 
                  "ep1_alzheimer_alg",                                
                  #"ep1_alzheimer_man", 
                  
                  "fup_vascular_dementia_yrs_alg",                
                  #"fup_vascular_dementia_yrs_man", 
                  "ep1_vascular_dementia_alg",                       
                  #"ep1_vascular_dementia_man", 
                  
                  "fup_stroke_yrs_alg",                               
                  #"fup_stroke_yrs_man", 
                  "ep1_stroke_alg",                                   
                  #"ep1_stroke_man", 
                  
                  "fup_ischstroke_yrs_alg",
                  #"fup_ischstroke_yrs_man", 
                  "ep1_ischstroke_alg"  )]                             
                  #"ep1_ischstroke_man")]

dim(df_use)
#502309     69

library(dplyr)
df_use = rename(df_use, eid = f.eid)

library(data.table)
bl_eid = as.data.frame(fread("UKB_NOV2023/Olink_BL_participants_50K_eid.txt.gz"))
df_olink_ppts = merge(df_use, bl_eid, by = "eid")

dim(df_olink_ppts)
#53022    69

fwrite(df_olink_ppts, 
       "UKB_NOV2023/Olink_BL_participants_50K_wOutcomes_wCog.txt.gz", 
       col.names=T, row.names=F, quote=F, sep="\t", append = F)

#=====================================================================================================================

# load outcome dataset 
df_olink_ppts = fread("Olink_BL_participants_50K_wOutcomes_wCog.txt.gz")

# load exposure dataset
ukb_rdata = readRDS("UK_Biobank/phenotype/UK_Biobank_Olink/olink_data_reformatted_081123.rds")
#dim(ukb_rdata)
# 55323  2925

# get the baseline measure 
ukb_bl = ukb_rdata[which(ukb_rdata$ins_index == 0), ]  # 53026
#dim(ukb_bl)
#53026  2925

# get the causal candidates 
ukb_bl_candidate = ukb_bl[,c("eid", "APOE", "MERTK", "TIMD4", "FLT4", "PDE5A", "HEXIM1", "CD46", 
                             "PEAR1", "EPHA2", "MEGF10", "NPTX1", "METAP1D", "TRIM21", "COL2A1")]

# merge exposure and outcome datasets
df_olink_exp_out = merge(df_olink_ppts, ukb_bl_candidate, by = "eid")
dim(df_olink_exp_out)
#53022    83

fwrite(df_olink_exp_out, 
       "UKB_NOV2023/Olink_BL_50K_wExpOutCog.txt.gz", 
       col.names=T, row.names=F, quote=F, sep="\t", append = F)

#===================================================================================================
# Further process the chosen baseline characteristics of this cohort 
#===================================================================================================

library(table1)
library(data.table)
library(dplyr)
#rm(list=ls())

df = fread("UKB_NOV2023/dataset/Olink_BL_50K_wExpOutCog.txt.gz")
df_diabetes_bl = fread("UKB_NOV2023/UKB_Diabetes1_2_baseline_22052024.txt")
df_diabetes_bl = df_diabetes_bl[ , -c("diabetes_t2d_selfreport", "date_baseline")]    # this variable is already included in the original df
  
# add in the baseline diabetic status as a new covariate
df = merge(df, df_diabetes_bl, by.x = "eid", by.y = "f.eid")
# all samples in the original df has a diabetic status (Yes/No) at baseline

#df = df_olink_exp_out

#ukb$edu_cat_detailed <- factor(ukb$edu_cat_detailed, levels=c("1", "2", "3", "4", "5", "6", "7"), 
#	labels=c("College or Univetabrsity degree", "A levels/AS levels or equivalent", "O levels/GCSEs or equivalent", "CSEs or equivalent", "NVQ or HND or HNC or equivalent", "Other professional qualifications (e.g. nursing, teaching)", "None of the above (primary education)"), ordered=FALSE)

#table(df$apoe_isoform_detailed) 
#table(df$apoe_isoform)
#table(df$apoe_e4_carrier)

df[,edu_yrs_imp := ifelse(edu_cat_detailed %in% c(7), 7,
                          ifelse(edu_cat_detailed %in% c(3, 4), 10,
                                 ifelse(edu_cat_detailed %in% c(2), 13,
                                        ifelse(edu_cat_detailed %in% c(5), 18,
                                               ifelse(edu_cat_detailed %in% c(1, 6), 20, as.integer(NA))))))]

df$sex_fac <- factor(df$sex, levels=c("0", "1"), labels=c("Female", "Male"), ordered=FALSE)
df$smallstat_fac <- factor(df$smallstat, levels=c("0", "1", "2"), labels=c("Never", "Previous", "Current"), ordered=FALSE)
df$alcallstat_fac <- factor(df$alcallstat, levels=c("0", "1", "2"), labels=c("Never", "Previous", "Current"), ordered=FALSE)

df$ep1_allcause_dementia_alg_fac <- factor(df$ep1_allcause_dementia_alg, levels=c("0", "1"), labels=c("Controls", "Incident dementia"), ordered=FALSE)
df$ep1_alzheimer_alg_fac <- factor(df$ep1_alzheimer_alg, levels=c("0", "1"), labels=c("Controls", "Incident AD"), ordered=FALSE)
df$ep1_vascular_dementia_alg_fac <- factor(df$ep1_vascular_dementia_alg, levels=c("0", "1"), labels=c("Controls", "Incident VaD"), ordered=FALSE)
df$ep1_stroke_alg_fac <- factor(df$ep1_stroke_alg, levels=c("0", "1"), labels=c("Controls", "Incident stroke"), ordered=FALSE)
df$ep1_ischstroke_alg_fac <- factor(df$ep1_ischstroke_alg, levels=c("0", "1"), labels=c("Controls", "Incident ischemic stroke"), ordered=FALSE)

#---------------------------------------
# 
# mannually check covariate distr.
#
#---------------------------------------

### categorical variables ###

table(df$sex_fac, useNA = "always")
table(df$smallstat_fac, useNA = "always")
table(df$alcallstat_fac, useNA = "always")
table(df$diabetes_t2d, useNA = "always")
table(df$apoe_e4_carrier, useNA = "always")

#table(df$sex_fac, useNA = "always")
#Female   Male   <NA> 
#28571  24420     31 

#table(df$smallstat_fac, useNA = "always")
#Never Previous  Current     <NA> 
#28680    18501     5603      238 

#table(df$alcallstat_fac, useNA = "always")
#Never Previous  Current     <NA> 
#2595     2197    48089      141 

#table(df$diabetes_t2d, useNA = "always")
#0     1  <NA> 
#47292  5729     1 

#table(df$apoe_e4_carrier, useNA = "always")
#0     1  <NA> 
#39771 13251     0 

table(df$ep1_allcause_dementia_alg, useNA = "always")
table(df$ep1_alzheimer_alg, useNA = "always")
table(df$ep1_vascular_dementia_alg, useNA = "always")
table(df$ep1_stroke_alg, useNA = "always")
table(df$ep1_ischstroke_alg, useNA = "always")

#table(df$ep1_allcause_dementia_alg, useNA = "always")
#0     1  <NA> 
#51551  1471     0 

#table(df$ep1_alzheimer_alg, useNA = "always")
#0     1  <NA> 
#52291   731     0 

#table(df$ep1_vascular_dementia_alg, useNA = "always")
#0     1  <NA> 
#52739   283     0 

#table(df$ep1_stroke_alg, useNA = "always")
#0     1  <NA> 
#51622  1400     0 

#table(df$ep1_ischstroke_alg, useNA = "always")
#0     1  <NA> 
#51829  1193     0 

### continuous variables ###

summary(df$age_baseline)
summary(df$edu_yrs_imp)
summary(df$townsend)
summary(df$bmi)
summary(df$sbp)
summary(df$tchol)
summary(df$hdl)
summary(df$ldl)
summary(df$crp)
summary(df$fram_points)
summary(df$duration1)

summary(df$fup_allcause_dementia_yrs_alg)
summary(df$fup_alzheimer_yrs_alg)
summary(df$fup_vascular_dementia_yrs_alg)
summary(df$fup_stroke_yrs_alg)
summary(df$fup_ischstroke_yrs_alg)

#summary(df$age_baseline)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#39.00   50.00   58.00   56.81   64.00   70.00 

#summary(df$edu_yrs_imp)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#7.0    10.0    18.0    14.7    20.0    20.0     838 

#summary(df$townsend)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-6.260  -3.620  -2.050  -1.179   0.770  10.420      66 

#summary(df$bmi)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#14.28   24.19   26.78   27.47   29.92   68.95     244 

#summary(df$sbp)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#76.5   124.5   136.5   137.8   149.5   246.5      44 

#summary(df$tchol)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#1.432   4.853   5.611   5.655   6.398  12.641    2462 

#summary(df$hdl)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.407   1.163   1.394   1.443   1.673   4.087    6594 

#summary(df$ldl)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.670   2.905   3.486   3.528   4.103   8.513    2568 

#summary(df$crp)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.080   0.660   1.350   2.662   2.820  79.960    2611 

#summary(df$fram_points)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-3.00    8.00   12.00   11.74   16.00   29.00     804 

#summary(df$duration1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#13.14   14.14   14.82   14.83   15.50   17.68 

#summary(df$fup_allcause_dementia_yrs_alg)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.03012 14.08624 14.75702 14.66711 15.48802 17.67556 

#summary(df$fup_alzheimer_yrs_alg)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.178  14.116  14.793  14.751  15.493  17.676 

#summary(df$fup_vascular_dementia_yrs_alg)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.345  14.133  14.806  14.801  15.500  17.676 

#summary(df$fup_stroke_yrs_alg)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.01643 14.08898 14.76249 14.65341 15.48802 17.67556 

#summary(df$fup_ischstroke_yrs_alg)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.01643 14.09446 14.76797 14.67931 15.48802 17.67556 

#summary(df_tab$cog_reaction_time)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-4.2117 -0.6451 -0.1546  0.0258  0.4855 11.3764     426 

#summary(df_tab$cog_numeric_memory)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-3.55   -0.52    0.24    0.02    0.99    4.02   38322 

#summary(df_tab$cog_fluid_intelligenc)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-2.830  -0.497  -0.031   0.015   0.902   3.235   26394 

#summary(df_tab$cog_pairs_matching)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-1.196043 -0.627931 -0.059819  0.001825  0.508293 19.540047 

#===============================================================
# Apply exclusion criteria
#===============================================================

# Step 1: Remove prevalent all-cause dementia and all-cause stroke cases
# Step 2: Remove ppts with NA in all chosen covariates 

df$dx_time_dementia = as.numeric(df$alg_allcause_dementia_date - df$date_baseline)
length(which(df$dx_time_dementia<=0))
# 48

df$dx_time_ad = as.numeric(df$alg_alzheimer_date - df$date_baseline)
length(which(df$dx_time_ad<=0))
# 6

df$dx_time_VaD = as.numeric(df$alg_vascular_dementia_date - df$date_baseline)
length(which(df$dx_time_VaD<=0))
# 1

df$dx_time_stroke = as.numeric(df$alg_stroke_date - df$date_baseline)
length(which(df$dx_time_stroke<=0))
# 902

df$dx_time_ischstroke = as.numeric(df$alg_ischstroke_date - df$date_baseline)
length(which(df$dx_time_ischstroke<=0))
# 244

df2 = df # copy a new dataframe in case the removal process goes wrong (so that the df will be untouched).

df2 = df %>% 
  filter(dx_time_stroke > 0 | is.na(dx_time_stroke)) %>%   # select incident any stroke cases or controls with no stroke events 
  filter(dx_time_dementia > 0 | is.na(dx_time_dementia))   # select incident all-cause dementia cases or controls with no dementia events 

which(df2$dx_time_stroke <=0)
which(df2$dx_time_ischstroke <=0)
which(df2$dx_time_dementia <=0)
which(df2$dx_time_ad <=0)
which(df2$dx_time_VaD <=0)
# checked: the prevalent stroke and dementia cases are removed. 

dim(df)
# 53022  108

dim(df2)
# 52076  108

#print(53022 - 52076)
# 946 excluded due to having dementia or stroke before baseline

df2_1 = df2 %>% 
  filter(!is.na(cog_reaction_time) & !is.na(cog_pairs_matching))

dim(df2_1)
#51332   108

#print(52076 - 51332)
# 744 were excluded due to no reaction time or pairs matching tests

df3 = df2_1 %>% 
  filter(!is.na(age_baseline) & !is.na(sex_fac) & !is.na(edu_yrs_imp) & !is.na(townsend)
         & !is.na(bmi) & !is.na(smallstat_fac) & !is.na(alcallstat_fac) & !is.na(sbp)
         & !is.na(tchol) & !is.na(ldl) & !is.na(crp) & !is.na(diabetes_bl)
         & !is.na(fram_points) & !is.na(apoe_e4_carrier))

# numeric memory and fluid intelligence were not used since they have too many missingness

dim(df3)
# 47571    108 (after further exclude those without reaction time)

df3 <-df3 %>%
  mutate(meanAge = mean(df3$age_baseline)) %>%
  mutate(stop_age_dementia = age_baseline + fup_allcause_dementia_yrs_alg) %>% 
  mutate(stop_age_ad = age_baseline + fup_alzheimer_yrs_alg) %>% 
  mutate(stop_age_VaD = age_baseline + fup_vascular_dementia_yrs_alg) %>% 
  mutate(stop_age_stroke = age_baseline + fup_stroke_yrs_alg) %>% 
  mutate(stop_age_ischstroke = age_baseline + fup_ischstroke_yrs_alg)

#print(52076 - 47571)  #4505
# 4505/52076 = 0.08650818 = 8.7%

# DONE with data processing #

fwrite(df3, 
       "UKB_NOV2023/dataset/Olink_BL_47571_wExpOutCog_14markers_cleaned_22052024.txt.gz", 
       col.names=T, row.names=F, quote=F, sep="\t", append = F)


################
# Table 1 ######
################

#--------------------------------
#
# further revision with Rstudio 
# add cog scores
#
#--------------------------------

library(table1)
library(data.table)

df_tab = fread("Olink_BL_47571_wExpOutCog_14markers_cleaned_22052024.txt.gz")
dim(df_tab)
# 47571   114

cog_outcome = fread("UKB_NOV2023/UKB_RT_PM_raw_baseline_27052024_v2.txt")
cog_out = cog_outcome %>% select(f.eid, rt_ms, pm_error_round2, rt_ln) %>% rename(eid = f.eid)

df_tab = inner_join(df_tab, cog_out)

which(is.na(df_tab$rt_ms))                   #none
which(is.na(df_tab$pm_error_round2))         #none

fwrite(df_tab, 
       "UKB_NOV2023/dataset/Olink_BL_47571_wExpOutCog_14markers_cleaned_27052024.txt.gz", 
       col.names=T, row.names=F, quote=F, sep="\t", append = F)

# cog_reaction_time and cog_pairs_matching in Eric's codes were Z-transformed; they did not follow normal distribution
# the raw value of the reaction time and the errors in pairs matching test was reported. 

df_tab = fread("Olink_BL_47571_wExpOutCog_14markers_cleaned_27052024.txt.gz")

# for smoking and alcohol drinking: 0-Never, 1-Previous, 2-Current
#df_tab$smallstat_bi = ifelse(df_tab$smallstat == 2, 1, 0)
#df_tab$alcallstat_bi = ifelse(df_tab$alcallstat == 2, 1, 0)

#df_tab$smallstat_bi_fac <- factor(df_tab$smallstat_bi, levels=c("0", "1"), labels=c("Not current", "Current"), ordered=FALSE)
#df_tab$alcallstat_bi_fac <- factor(df_tab$alcallstat_bi, levels=c("0", "1"), labels=c("Not current", "Current"), ordered=FALSE)

# all-cause dementia
table1(~ age_baseline + sex_fac + edu_yrs_imp + townsend + bmi + smallstat_bi_fac + alcallstat_bi_fac + sbp + 
             tchol +  ldl + diabetes_bl_fac + as.factor(apoe_e4_carrier) + rt_ms + pm_error_round2 +
             duration1 + fup_allcause_dementia_yrs_alg | ep1_allcause_dementia_alg_fac, dat = df_tab)                  

summary(df_tab[which(df_tab$ep1_allcause_dementia_alg_fac == "Controls"), ]$rt_ms)
summary(df_tab[which(df_tab$ep1_allcause_dementia_alg_fac == "Incident dementia"), ]$rt_ms)

summary(df_tab[which(df_tab$ep1_allcause_dementia_alg_fac == "Controls"), ]$pm_error_round2)
summary(df_tab[which(df_tab$ep1_allcause_dementia_alg_fac == "Incident dementia"), ]$pm_error_round2)


# all stroke
table1(~ age_baseline + sex_fac + edu_yrs_imp + townsend + bmi + smallstat_bi_fac + alcallstat_bi_fac + sbp + 
            tchol +  ldl + diabetes_bl_fac + as.factor(apoe_e4_carrier) + rt_ms + pm_error_round2 + 
            duration1 + fup_stroke_yrs_alg | ep1_stroke_alg_fac, dat = df_tab)  

summary(df_tab[which(df_tab$ep1_stroke_alg_fac == "Controls"), ]$rt_ms)
summary(df_tab[which(df_tab$ep1_stroke_alg_fac == "Incident stroke"), ]$rt_ms)

summary(df_tab[which(df_tab$ep1_stroke_alg_fac == "Controls"), ]$pm_error_round2)
summary(df_tab[which(df_tab$ep1_stroke_alg_fac == "Incident stroke"), ]$pm_error_round2)

# total pop 
summary(df_tab$rt_ms)
summary(df_tab$pm_error_round2)

# get the years of FU --> incident dementia, stroke, or censoring whichever occurs first
df_tab$overall_FU = pmin(df_tab$fup_allcause_dementia_yrs_alg, df_tab$fup_stroke_yrs_alg, na.rm = T)

hist(df_tab$age_baseline, 
     xlab = "Age at baseline assessment", 
     ylab = "Number of participants",
     main = "Histogram of age distribution at UKB baseline")

hist(df_tab$rt_ms, 
     xlab = "Reaction time at baseline (ms)", 
     ylab = "Number of participants")

hist(df_tab$rt_ln, 
     xlab = "log reaction time at baseline", 
     ylab = "Number of participants")

hist(df_tab$pm_error_round2, 
     xlab = "Errors in pairs matching test Round 2 at baseline", 
     ylab = "Number of participants")

df_tab$pm_error_round2[which(df_tab$pm_error_round2 > 30)]

#########################################################################
# based on the z-transformed result, outliers still exist for some markers (see plots)
# consider to use inverse normal transformation 
#########################################################################

df_tab = df_tab %>%
  mutate(APOE_int = qnorm((rank(APOE,na.last="keep")-0.5)/sum(!is.na(APOE)))) %>%
  mutate(MERTK_int = qnorm((rank(MERTK,na.last="keep")-0.5)/sum(!is.na(MERTK)))) %>%
  mutate(TIMD4_int = qnorm((rank(TIMD4,na.last="keep")-0.5)/sum(!is.na(TIMD4)))) %>%
  mutate(FLT4_int = qnorm((rank(FLT4,na.last="keep")-0.5)/sum(!is.na(FLT4)))) %>%
  mutate(PDE5A_int = qnorm((rank(PDE5A,na.last="keep")-0.5)/sum(!is.na(PDE5A)))) %>%
  mutate(HEXIM1_int = qnorm((rank(HEXIM1,na.last="keep")-0.5)/sum(!is.na(HEXIM1)))) %>%
  mutate(CD46_int = qnorm((rank(CD46,na.last="keep")-0.5)/sum(!is.na(CD46)))) %>%
  mutate(PEAR1_int = qnorm((rank(PEAR1,na.last="keep")-0.5)/sum(!is.na(PEAR1)))) %>%
  mutate(EPHA2_int = qnorm((rank(EPHA2,na.last="keep")-0.5)/sum(!is.na(EPHA2)))) %>%
  mutate(MEGF10_int = qnorm((rank(MEGF10,na.last="keep")-0.5)/sum(!is.na(MEGF10)))) %>%
  mutate(NPTX1_int = qnorm((rank(NPTX1,na.last="keep")-0.5)/sum(!is.na(NPTX1)))) %>%
  mutate(METAP1D_int = qnorm((rank(METAP1D,na.last="keep")-0.5)/sum(!is.na(METAP1D)))) %>%
  mutate(TRIM21_int = qnorm((rank(TRIM21,na.last="keep")-0.5)/sum(!is.na(TRIM21)))) %>%
  mutate(COL2A1_int = qnorm((rank(COL2A1,na.last="keep")-0.5)/sum(!is.na(COL2A1))))
  
## RE-PLOT histograms
df_tab = as.data.frame(df_tab)

pdf("UKB_NOV2023/Cog_Inflam_INT_distribution_UKBPPP_BL_cohort_remove_OutCovMissing.pdf", width=10,height=5)
par(mfrow = c(1, 2))

for(i in c(118:131)) {
  print(hist(df_tab[ ,i], main = "", xlab = colnames(df_tab)[i]))
}
dev.off()         

###################################
# save data 

#fwrite(df_tab, 
#       "UKB_NOV2023/Olink_BL_47571_wExpOutCog_INT_14markers_28052024.txt.gz", 
#       col.names=T, row.names=F, quote=F, sep="\t", append = F)

df_tab = fread("UKB_NOV2023/Olink_BL_47571_wExpOutCog_INT_14markers_28052024.txt.gz") 

# 0-Never, 1-Previous, 2-Current 
# binary categorization: 0-Not current, 1-Current
df_tab$smallstat_bi = ifelse(df_tab$smallstat == 2, 1, 0)
df_tab$alcallstat_bi = ifelse(df_tab$alcallstat == 2, 1, 0)

df_tab$smallstat_bi_fac <- factor(df_tab$smallstat_bi, levels=c("0", "1"), labels=c("Not current", "Current"), ordered=FALSE)
df_tab$alcallstat_bi_fac <- factor(df_tab$alcallstat_bi, levels=c("0", "1"), labels=c("Not current", "Current"), ordered=FALSE)

fwrite(df_tab, 
       "UKB_NOV2023/Olink_BL_47571_wExpOutCog_INT_14markers_06062024.txt.gz", 
       col.names=T, row.names=F, quote=F, sep="\t", append = F)

##################################################################################################################################################

#---------------------
#
# survival analysis
#
#---------------------

# Perform analysis for all-cause dementia and stroke 

library(data.table)
library(dplyr)
library(survival)
#rm(list=ls())

df = fread("UKB_NOV2023/Olink_BL_47571_wExpOutCog_INT_14markers_06062024.txt.gz")
dim(df)
#47571  135

#==========================================================================================================================================
#
# Fit Cox proportional hazards models
# 
# Age at baseline recruitment is used as the underlying timescale. 
# 
# For all-cause dementia 
# Model1: adjust for sex, education in years, APOE e4 carrier status
# Model2: adjust for sex, education in years, APOE e4 status, Townsend deprevation index (continuous), 
#         BMI (continuous), smoking status (past, previous, current), alcohol drinking (past, previous, current), 
#         SBP, total cholesterol, LDL, baseline diabetes (Type I and II). 
#
#==========================================================================================================================================

#--------------------------
#
# for dementia modeling 
#
#--------------------------

# Model 1: age as the timescale, adjusting for sex, education, and apoe e4 status

biomarker_list = colnames(df)[c(118:131)]

# i=1 for test
for(i in 2:length(biomarker_list)) {
  
  tryCatch({
    
    age.m1.acd <-coxph(Surv(age_baseline,stop_age_dementia, ep1_allcause_dementia_alg==1) ~ 
                         get(biomarker_list[i]) + sex + edu_yrs_imp + apoe_e4_carrier, data = df, ties='efron')
    
    res.m1.acd = summary(age.m1.acd)
    saveRDS(res.m1.acd, paste0("UKB_NOV2023/proteomic_14_cox_reg/", biomarker_list[i] , "_m1_res_acd_06062024.RDS"))
    
    res_ln = as.data.frame(res.m1.acd$coefficients)[1, ]
    conf_int = as.data.frame(res.m1.acd$conf.int)[1, ]
    
    res_line = merge(res_ln, conf_int)
    res_line$n = res.m1.acd$n
    res_line$n_event = res.m1.acd$nevent
    res_line$exposure = biomarker_list[i]
    res_line$outcome = "all_cause_dementia"
    res_line$model = "primary"
    
    fwrite(res_line, "UKB_NOV2023/proteomic_14_cox_reg/cox_14markers_res_06062024.txt", col.names=F, row.names=F, quote=F, sep="\t", append = T)

      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# Model 2: age as the timescale, adjusting for sex, education, and apoe e4 status, additionally adjust for other vascular related risk factors 

for(i in 1:length(biomarker_list)) {
  
  tryCatch({
    
    age.m2.acd <-coxph(Surv(age_baseline,stop_age_dementia, ep1_allcause_dementia_alg==1) ~ 
                         get(biomarker_list[i]) + sex + edu_yrs_imp + apoe_e4_carrier + 
                         townsend + bmi + smallstat_bi + alcallstat_bi + sbp + 
                         tchol +  ldl + diabetes_bl, data = df, ties='efron')
    
    res.m2.acd = summary(age.m2.acd)
    saveRDS(res.m2.acd, paste0("UKB_NOV2023/proteomic_14_cox_reg/", biomarker_list[i] , "_m2_res_acd_06062024.RDS"))
    
    res_ln = as.data.frame(res.m2.acd$coefficients)[1, ]
    conf_int = as.data.frame(res.m2.acd$conf.int)[1, ]
    
    res_line = merge(res_ln, conf_int)
    res_line$n = res.m2.acd$n
    res_line$n_event = res.m2.acd$nevent
    res_line$exposure = biomarker_list[i]
    res_line$outcome = "all_cause_dementia"
    res_line$model = "secondary"
    
    fwrite(res_line, "UKB_NOV2023/proteomic_14_cox_reg/cox_14markers_res_06062024.txt", col.names=F, row.names=F, quote=F, sep="\t", append = T)
  
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#--------------------------
#
# for any stroke modeling 
#
#--------------------------

# Model 1: age as the timescale, adjusting for sex

for(i in 1:length(biomarker_list)) {
  
  tryCatch({
    age.m1.stroke <-coxph(Surv(age_baseline,stop_age_stroke, ep1_stroke_alg==1) ~ 
                            get(biomarker_list[i]) + sex, data = df, ties='efron')
    
    res.m1.stroke = summary(age.m1.stroke)
    saveRDS(res.m1.stroke, paste0("UKB_NOV2023/proteomic_14_cox_reg/", biomarker_list[i] , "_m1_res_as_06062024.RDS"))
    
    res_ln = as.data.frame(res.m1.stroke$coefficients)[1, ]
    conf_int = as.data.frame(res.m1.stroke$conf.int)[1, ]
    
    res_line = merge(res_ln, conf_int)
    res_line$n = res.m1.stroke$n
    res_line$n_event = res.m1.stroke$nevent
    res_line$exposure = biomarker_list[i]
    res_line$outcome = "any_stroke"
    res_line$model = "primary"
    
    fwrite(res_line, "UKB_NOV2023/proteomic_14_cox_reg/cox_14markers_res_06062024.txt", col.names=F, row.names=F, quote=F, sep="\t", append = T)
  
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# Model 2: age as the timescale, adjusting for sex, additionally adjust for other vascular related risk factors 

for(i in 1:length(biomarker_list)) {
  
  tryCatch({
    
    age.m2.stroke <-coxph(Surv(age_baseline,stop_age_stroke, ep1_stroke_alg==1) ~ 
                            get(biomarker_list[i]) + sex +
                            townsend + bmi + smallstat_bi + alcallstat_bi + sbp + 
                            tchol +  ldl + diabetes_bl, data = df, ties='efron')
    
    res.m2.stroke = summary(age.m2.stroke)
    saveRDS(res.m2.stroke, paste0("UKB_NOV2023/proteomic_14_cox_reg/", biomarker_list[i] , "_m2_res_as_06062024.RDS"))
    
    res_ln = as.data.frame(res.m2.stroke$coefficients)[1, ]
    conf_int = as.data.frame(res.m2.stroke$conf.int)[1, ]
    
    res_line = merge(res_ln, conf_int)
    res_line$n = res.m2.stroke$n
    res_line$n_event = res.m2.stroke$nevent
    res_line$exposure = biomarker_list[i]
    res_line$outcome = "any_stroke"
    res_line$model = "secondary"
    
    fwrite(res_line, "UKB_NOV2023/proteomic_14_cox_reg/cox_14markers_res_06062024.txt", col.names=F, row.names=F, quote=F, sep="\t", append = T)
  
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


#==========================================================================================================================================
#
# Fit linear regression models
# 
# For all-cause dementia 
# Model1: adjust for age (linear + quadratic), sex, education in years, APOE e4 carrier status
# Model2: adjust for age (linear + quadratic), sex, education in years, APOE e4 status, Townsend deprivation index (continuous),
#         BMI (continuous), smoking status (past, previous, current), alcohol drinking (past, previous, current), 
#         SBP, total cholesterol, LDL, Type II diabetes. 
#
#==========================================================================================================================================

#########################################
# Instance 0 -- Baseline 
# outcomes: reaction time at baseline
#########################################

#------------
#
# Model basic
#
#------------

for(i in 2:length(biomarker_list)) {
  tryCatch({
    
    m.rt.bl <- glm(rt_ln ~ get(biomarker_list[i]) + sex + age_baseline + I(age_baseline^2), data = df)
    
    summ = summary(m.rt.bl)
    
    beta <- coef(m.rt.bl)[2]
    SE <- coef(summary(m.rt.bl))[2,2]
    lcl <- beta-qnorm(0.975)*SE 
    ucl <- beta+qnorm(0.975)*SE
    t_stats = coef(summary(m.rt.bl))[2,3]
    p_val = coef(summary(m.rt.bl))[2,4]
    
    summ_tab = as.data.frame(cbind(beta, lcl, ucl, p_val, SE, t_stats))
    summ_tab$n = summ$df[2] +5                       # n of obs = df + 5 variables (incl. 4 independent variables + 1 intercept)
    summ_tab$exposure = biomarker_list[i]
    summ_tab$outcome = "natural log of reaction time"
    summ_tab$model = "primary"
    
    #tot_tab = summ_tab
    tot_tab = rbind(tot_tab, summ_tab)                      # run the first marker in the loop to get first tot_tab (i.e. tot_tab = summ_tab for the 1st marker) then loop the rest
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.table(tot_tab, 
            "UKB_NOV2023/proteomic_14_linear_reg/reaction_time_14markers_res_06062024.txt", col.names=T, row.names=F, quote=F, sep="\t", append = F)

#------------
#
# Model VRF
#
#------------

for(i in 2:length(biomarker_list)) {
  tryCatch({
    
    m.rt.bl <- glm(rt_ln ~ get(biomarker_list[i]) + sex + age_baseline + I(age_baseline^2) + edu_yrs_imp + I(edu_yrs_imp^2) + 
                     townsend + I(townsend^2) + bmi + I(bmi^2) + smallstat_bi + alcallstat_bi + sbp + 
                     tchol + ldl + diabetes_bl , data = df)
    
    summ = summary(m.rt.bl)

    beta <- coef(m.rt.bl)[2]
    SE <- coef(summary(m.rt.bl))[2,2]
    lcl <- beta-qnorm(0.975)*SE 
    ucl <- beta+qnorm(0.975)*SE
    t_stats = coef(summary(m.rt.bl))[2,3]
    p_val = coef(summary(m.rt.bl))[2,4]
    
    summ_tab = as.data.frame(cbind(beta, lcl, ucl, p_val, SE, t_stats))
    summ_tab$n = summ$df[2] + 17 
    summ_tab$exposure = biomarker_list[i]
    summ_tab$outcome = "natural log of reaction time"
    summ_tab$model = "secondary"
    
    #tot_tab = summ_tab
    tot_tab = rbind(tot_tab, summ_tab)                      # run the first marker in the loop to get first tot_tab (i.e. tot_tab = summ_tab for the 1st marker) then loop the rest
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.table(tot_tab, 
            "UKB_NOV2023/proteomic_14_linear_reg/reaction_time_14markers_res_06062024.txt", col.names=F, row.names=F, quote=F, sep="\t", append = T)


##############################################
# Instance 0 -- Baseline 
# outcomes: cog_pairs_matching at baseline
##############################################

# try negative binomial modeling vs. ln(score+1) modeling

# because many ppts had no errors in the test, so we have to plus 1 to do the natural log transformation 
df_tab$pm_ln = log(df_tab$pm_error_round2 + 1)

hist(df_tab$pm_ln, 
     xlab = "log pairs matching errors at baseline", 
     ylab = "Number of participants")

# the histogram does not look "normal" after the natural log transformation and the events are discrete 
# thus, we decide to try the negative binomial methods. 

#------------
#
# Model basic
#
#------------

#install.packages("MASS")
library(MASS)

for(i in 2:length(biomarker_list)) {
  tryCatch({
    
    m.rt.bl <- glm.nb(pm_error_round2 ~ get(biomarker_list[i]) + sex + age_baseline, data = df)
    
    summ = summary(m.rt.bl)
    
    beta <- coef(m.rt.bl)[2]
    SE <- coef(summary(m.rt.bl))[2,2]
    lcl <- confint(m.rt.bl)[2,1]
    ucl <- confint(m.rt.bl)[2,2]
    p_val = coef(summary(m.rt.bl))[2,4]
    RR = exp(beta)
    RR_lcl = exp(lcl)
    RR_ucl = exp(ucl)
    
    summ_tab = as.data.frame(cbind(beta, lcl, ucl, p_val, SE, RR, RR_lcl, RR_ucl))
    summ_tab$n = summ$df[2]+4    # n of obs = df + 3 independent variables + 1 intercept/dispersion parameter 
    summ_tab$exposure = biomarker_list[i]
    summ_tab$outcome = "errors in pairs matching"
    summ_tab$model = "primary"
    
    #tot_tab = summ_tab
    tot_tab = rbind(tot_tab, summ_tab)                      # run the first marker in the loop to get first tot_tab (i.e. tot_tab = summ_tab for the 1st marker) then loop the rest
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.table(tot_tab, 
            "UKB_NOV2023/proteomic_14_linear_reg/pairs_matching_14markers_res_06062024.txt", col.names=T, row.names=F, quote=F, sep="\t", append = F)

#------------
#
# Model VRF
#
#------------

for(i in 2:length(biomarker_list)) {
  tryCatch({
    
   m.rt.bl <- glm.nb(pm_error_round2 ~ get(biomarker_list[i]) + sex + age_baseline + edu_yrs_imp + 
                        townsend + bmi + smallstat_bi + alcallstat_bi + sbp + 
                        tchol + ldl + diabetes_bl, data = df)
    
    summ = summary(m.rt.bl)
    
    beta <- coef(m.rt.bl)[2]
    SE <- coef(summary(m.rt.bl))[2,2]
    lcl <- confint(m.rt.bl)[2,1]
    ucl <- confint(m.rt.bl)[2,2]
    p_val = coef(summary(m.rt.bl))[2,4]
    RR = exp(beta)
    RR_lcl = exp(lcl)
    RR_ucl = exp(ucl)
    
    summ_tab = as.data.frame(cbind(beta, lcl, ucl, p_val, SE, RR, RR_lcl, RR_ucl))
    summ_tab$n = summ$df[2]+13
    summ_tab$exposure = biomarker_list[i]
    summ_tab$outcome = "errors in pairs matching"
    summ_tab$model = "secondary"
    
    #tot_tab = summ_tab
    tot_tab = rbind(tot_tab, summ_tab)                     # run the first marker in the loop to get first tot_tab (i.e. tot_tab = summ_tab for the 1st marker) then loop the rest
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.table(tot_tab, 
            "UKB_NOV2023/proteomic_14_linear_reg/pairs_matching_14markers_res_06062024.txt", col.names=F, row.names=F, quote=F, sep="\t", append = T)


#--------------------
#
# Result organization 
#
#--------------------

# based on result examination, the signals are consistent across NPX values of exposure, Z-transformed NPX values, and inverse rank transformed NPX values 
# in addition, the results from basically adjusted model and VRF adjusted model are largely consistent 
# for illustration purpose, we used results from basically adjusted model (M1) and the INT transformed values 

#----------------
#
# forest plot 
#
#----------------

setwd("Supplemental Materials/original_data_sheet")

library(data.table)
library(dplyr)
#library(stringr)
library(tidyr) 

res_cox = fread("cox_14_res_org.txt")
# the results are consist across raw NPX values, Z-transformed NPX, and INT NPX as well as across basic and fully adjusted models. 
# basic model (m1) and inverse normal transformed version is used as the primary results

res_cox = res_cox %>% 
  filter(exposure != "TRIM21_int") %>% 
  mutate(report_est = paste0(format(round(HR, digits = 2), nsmall=2), 
                             " (", format(round(lcl, digits = 2), nsmall=2), ", ", format(round(ucl, digits = 2), nsmall=2), ") ")) 

acd_m1_cox = res_cox %>% 
  filter(outcome == "all_cause_dementia" & model == "primary") %>%
  mutate(p_fdr = p.adjust(p_val, method = "BH"))
  
colnames(acd_m1_cox) = c("model", "exposure", "outcome.acd", "beta.acd", "se.acd", "p_val.acd",
                         "HR.acd", "lcl.acd", "ucl.acd", "report_est.acd", "p_fdr.acd") 

as_m1_cox = res_cox %>% 
  filter(outcome == "any_stroke" & model == "primary") %>%
  mutate(p_fdr = p.adjust(p_val, method = "BH"))

colnames(as_m1_cox) = c("model", "exposure", "outcome.as", "beta.as", "se.as", "p_val.as",
                         "HR.as", "lcl.as", "ucl.as", "report_est.as", "p_fdr.as") 


df_wide = inner_join(acd_m1_cox, as_m1_cox)
df_wide = df_wide %>% 
  separate(exposure, "Gene", sep = "_")

df_wide$label_name <- factor(df_wide$Gene, levels = c("HEXIM1",
                                                      "CD46",
                                                      "PEAR1",
                                                      "PDE5A",
                                                      "APOE", 
                                                      "NPTX1", 
                                                      "COL2A1",
                                                      "MERTK",
                                                      "FLT4",
                                                      "TIMD4",
                                                      "MEGF10",
                                                      "EPHA2", 
                                                      "METAP1D")) 

library(forestploter)
options(digits=3)

# Set-up theme
tm <- forest_theme(base_size = 11,
                   refline_col = "#222222",
                   legend_name = "Outcomes",
                   legend_value = c("All-cause dementia", "Any stroke"))


# plot (combined for dementia and stroke)
dt = df_wide[order(df_wide$label_name, decreasing = TRUE), ]
dt$space1 = paste(rep("    ", 28), collapse = " ")

p <- forest(dt[ ,c("label_name", "space1",  )],
            
            est = list(dt$HR.acd, 
                       dt$HR.as),
            
            lower = list(dt$lcl.acd, 
                         dt$lcl.as), 
            
            upper = list(dt$ucl.acd, 
                         dt$ucl.as),
            
            ci_column = 2,
            ref_line = 1,
            nudge_y = 0.3,
            xlim = c(0.2, 1.8),
            xlab = "HR",
            theme = tm,
            ticks_at = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8))

#plot(p)

#--------------------------
# forest plot for dementia 
#--------------------------

# Set-up theme
tm <- forest_theme(base_size = 11,
                   refline_col = "#222222", 
                   ci_col = "#BF40BF",  #
                   ci_fill = "#BF40BF")

# plot (dementia)
dt = df_wide[order(df_wide$label_name, decreasing = TRUE), ]
dt$space1 = paste(rep("    ", 28), collapse = " ")

p <- forest(dt[ ,c("label_name", "space1", "report_est.acd", "p_fdr.acd")],
            
            est = dt$HR.acd, 
            
            lower = dt$lcl.acd, 
            
            upper = dt$ucl.acd, 
            
            ci_column = 2,
            ref_line = 1,
            nudge_y = 0.3,
            xlim = c(0.5, 1.3),
            xlab = "HR",
            theme = tm,
            ticks_at = c(0.6, 0.8, 1.0, 1.2))

ggplot2::ggsave(filename = "forest_for_cox_dementia_28052024.pdf", plot = p,
                dpi = 600,
                width = 9.5, height = 4, units = "in")

#--------------------------
# forest plot for stroke 
#--------------------------

# Set-up theme
tm <- forest_theme(base_size = 11,
                   refline_col = "#222222", 
                   ci_col = "#f27a74",  ##BF40BF
                   ci_fill = "#f27a74")

# plot stroke
dt = df_wide[order(df_wide$label_name, decreasing = TRUE), ]
dt$space1 = paste(rep("    ", 28), collapse = " ")

p <- forest(dt[ ,c("label_name", "space1", "report_est.as", "p_fdr.as")],
            
            est = dt$HR.as, 
            
            lower = dt$lcl.as, 
            
            upper = dt$ucl.as, 
            
            ci_column = 2,
            ref_line = 1,
            nudge_y = 0.3,
            xlim = c(0.9, 1.31),
            xlab = "HR",
            theme = tm,
            ticks_at = c(0.9, 1.0, 1.2))

ggplot2::ggsave(filename = "forest_for_cox_stroke_28052024.pdf", plot = p,
                dpi = 600,
                width = 9.5, height = 4, units = "in")

#-------------------------------------------------------
#
# Check the 13 markers identified by MR 
# 
# Compare cross-sectional analysis on reaction time
#
#-------------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)      # the library used for separate function

res_cog_bl = fread("cog_perform_res_org.txt")

res_rt_bl = res_cog_bl[1:14, ]
res_pm_bl = res_cog_bl[29:42, ]

res_rt_bl = res_rt_bl %>% separate(exposure, "Gene", sep = "_") %>% filter(Gene != "TRIM21") %>% mutate(p_fdr = p.adjust(p_val, method = "BH"))
res_pm_bl = res_pm_bl %>% separate(exposure, "Gene", sep = "_") %>% filter(Gene != "TRIM21") %>% mutate(p_fdr = p.adjust(p_val, method = "BH"))

res_rt_bl = res_rt_bl %>% mutate(report_est = paste0(format(round(percent_change, digits = 2), nsmall=2), 
                           " (", format(round(percent_change_lcl, digits = 2), nsmall=2), ", ", format(round(percent_change_ucl, digits = 2), nsmall=2), ") ")) 

res_pm_bl = res_pm_bl %>% mutate(report_est = paste0(format(round(percent_change, digits = 2), nsmall=2), 
                                                     " (", format(round(percent_change_lcl, digits = 2), nsmall=2), ", ", format(round(percent_change_ucl, digits = 2), nsmall=2), ") ")) 


res_rt_bl$label_name <- factor(res_rt_bl$Gene, levels = c("HEXIM1",
                                                          "CD46",
                                                          "PEAR1",
                                                          "PDE5A",
                                                          "APOE", 
                                                          "NPTX1", 
                                                          "COL2A1",
                                                          "MERTK",
                                                          "FLT4",
                                                          "TIMD4",
                                                          "MEGF10",
                                                          "EPHA2", 
                                                          "METAP1D")) 


res_pm_bl$label_name <- factor(res_pm_bl$Gene, levels = c("HEXIM1",
                                                          "CD46",
                                                          "PEAR1",
                                                          "PDE5A",
                                                          "APOE", 
                                                          "NPTX1", 
                                                          "COL2A1",
                                                          "MERTK",
                                                          "FLT4",
                                                          "TIMD4",
                                                          "MEGF10",
                                                          "EPHA2", 
                                                          "METAP1D")) 

library(forestploter)
options(digits=3)

#-----------
#
# RT
#
#-----------

# Set-up theme
tm <- forest_theme(base_size = 11,
                   refline_col = "#222222",
                   ci_col = "#4169E1",  ##BF40BF
                   ci_fill = "#4169E1")

# plot
dt = res_rt_bl[order(res_rt_bl$label_name, decreasing = TRUE), ]
dt$space1 = paste(rep("    ", 28), collapse = " ")

p <- forest(dt[ ,c("label_name", "space1", "report_est", "p_fdr")],
            
            est = list(dt$percent_change),
            
            lower = list(dt$percent_change_lcl), 
            
            upper = list(dt$percent_change_ucl),
            
            ci_column = 2,
            ref_line = 0,
            nudge_y = 0.3,
            xlim = c(-0.8, 0.8),
            xlab = "Percent difference in reaction time",
            theme = tm,
            ticks_at = c(-0.7, -0.3, 0, 0.3, 0.7))
#plot(p)

ggplot2::ggsave(filename = "forest_for_reactionTime_28052024.pdf", plot = p,
                dpi = 600,
                width = 9.5, height = 4, units = "in")


#---------------
#
# Pairs matching
#
#---------------

# Set-up theme
tm <- forest_theme(base_size = 11,
                   refline_col = "#222222",
                   ci_col = "#008080",  ##BF40BF
                   ci_fill = "#008080")

# plot
dt = res_pm_bl[order(res_pm_bl$label_name, decreasing = TRUE), ]
dt$space1 = paste(rep("    ", 28), collapse = " ")
dt$p_fdr = format(dt$p_fdr, scientific = T)

p <- forest(dt[ ,c("label_name", "space1", "report_est", "p_fdr")],
            
            est = list(dt$percent_change),
            
            lower = list(dt$percent_change_lcl), 
            
            upper = list(dt$percent_change_ucl),
            
            ci_column = 2,
            ref_line = 0,
            nudge_y = 0.3,
            xlim = c(-2.3, 1.6),
            xlab = "Percent difference in the number of matching errors",
            theme = tm,
            ticks_at = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5))
#plot(p)

ggplot2::ggsave(filename = "forest_for_PairsMatch_28052024.pdf", plot = p,
                dpi = 600,
                width = 9.5, height = 4, units = "in")


