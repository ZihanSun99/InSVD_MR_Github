
#-----------------------------------------------
#
# Coloc between exposure and outcome GWAS 
#
# binary imaging outcomes (CMB, PVS, LS)
# continuous imaging outcomes (WMH, FA, MD)
#
# primary analysis: p12 = 1e-5
#-----------------------------------------------

library(data.table)
library(stringr)
library(dplyr)
library(coloc)
library(tidyr)

#rm(list=ls())

writeLines("Start Job Num")

job_id <- commandArgs(trailingOnly=TRUE)[1]
print(job_id)
job_id <- as.numeric(job_id)

#job_id = 1

#=====================================================================================================================================================
# make the array list

prot_tss = fread("ukb_prot_tss_949_04032024.txt") 

prot_tss = prot_tss %>% 
  mutate(Gene_name = folder_name) %>%
  separate(Gene_name, "Gene", sep = "_" )

prot_tss = prot_tss %>% 
  filter(Gene == "HEXIM1"|
           Gene =="CD46"|
           Gene =="PEAR1"|
           Gene =="PDE5A"|
           Gene =="APOE"| 
           Gene =="NPTX1"| 
           Gene =="COL2A1"|
           Gene == "MERTK"|
           Gene =="FLT4"|
           Gene =="TIMD4"|
           Gene =="MEGF10"|
           Gene == "EPHA2"| 
           Gene =="METAP1D")

outcomes = c("EPVS", "CMB")
#outcomes = "LS_TOAST_MRI"
#outcomes = c("WMH", "FA", "MD")

array_list  = prot_tss %>% 
  slice(rep(1:nrow(prot_tss), each = 2)) %>%     
  mutate(outcome = rep(outcomes, 13)) %>%   
  mutate(CHR = as.numeric(CHR)) %>%
  mutate(gene_start = as.numeric(gene_start)) %>%
  mutate(gene_end = as.numeric(gene_end))

#=====================================================================================================================================================

exposure = array_list$folder_name[job_id]
outcomes = array_list$outcome[job_id]
cis_chr = array_list$CHR[job_id]

lbound = array_list$gene_start[job_id] - 200000
lbound = ifelse(lbound < 0, 0, lbound)

ubound = array_list$gene_end[job_id] + 200000

print(paste0(cis_chr, ":", lbound, "-", ubound))

inflam_file = fread(paste0("coloc/exposure_cis_aligned_file/", exposure, "_cis_for_coloc.txt.gz"))
outcome_file  = fread(paste0("GWAS_Sumstats_Outcomes/hg38_version/by_chromosome/", outcomes, "_", cis_chr, "_GWAS_Sumstats_GRCh38_cleaned.txt.gz"))

inflam_200kb = inflam_file %>% 
  filter(BP >= lbound & BP <= ubound) 

common_snps = inflam_200kb %>% 
  select(VarID) %>% 
  inner_join(outcome_file) %>%
  select(VarID)

df_exp = inflam_200kb %>% 
  inner_join(common_snps) %>%
  mutate(VARBETA = SE^2)

df_out = outcome_file %>% 
  inner_join(common_snps) %>%
  mutate(VARBETA = SE^2)

#identical(df_exp$VarID, df_out$VarID)
#identical(df_exp$A0, df_out$A0)
#identical(df_exp$A1, df_out$A1)

#which(is.na(df_out$BETA))
#min(df_out$VARBETA)
#min(df_exp$VARBETA)

#==================================================================================================================================================================================================================
#minimum_data=D1[c("beta","varbeta","snp","position","type","sdY")]
#list=as.list(df)

df_exp_list=list(
  snp= df_exp$VarID, 
  beta= df_exp$BETA, 
  varbeta= df_exp$VARBETA, 
  position= df_exp$BP, 
  type= "quant", 
  sdY= 1)

df_out_list=list(
  snp= df_out$VarID, 
  beta= df_out$BETA, 
  varbeta= df_out$VARBETA, 
  position= df_out$BP, 
  type= "cc",              # if the outcome is WMH, FA, or MD, then change this to "quant"
  sdY= 1)

res = coloc.abf(df_exp_list, df_out_list, MAF = NULL, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)

# extract causal snps conditional on H4 is true
causal_snp = subset(res$results,SNP.PP.H4>0.01) 

# extract the 95% credible set
o <- order(res$results$SNP.PP.H4,decreasing=TRUE)
# step 1: order the SNP.PP.H4 in an descending order

cs <- cumsum(res$results$SNP.PP.H4[o])
# step 2: add the SNP.PP.H4 from its largest to smallest value
# cumsum returns a vector whose elements are the cumulative sums
# of the elements of the argument.

w <- which(cs > 0.95)[1]
# step 3: get the index of the first SNP whose PP.H4 exceeds 95

cs_snp = as.data.frame(res$results[o,][1:w,]$snp)
colnames(cs_snp) = "credible_set"
# step 4: obtain the credible set from the ordered SNP list up till this index SNP

#----------------------
# save all the results 
#----------------------

# save the full colocalization results: SNP Priors, Hypothesis Priors, Posterior
saveRDS(res, paste0("coloc/200kb/pp_res/", exposure, "_", outcomes, "_200kb_res_cc_v2.RDS"))

# save the posterior probability of the 4 Hypothesis 
res_sum = as.data.frame(t(res$summary))
res_sum$gwas.exposure = exposure
res_sum$gwas.outcome = outcomes
res_sum$PP.H4.threshold = 0.8
res_sum$colocalized = ifelse(res_sum$PP.H4.abf >= 0.8, "yes", "no")

fwrite(res_sum, "coloc/pp_h4_200kb_res_pvs_cmb_cc_v2_22052024.txt", col.names=F, row.names=F, quote=F, sep="\t", append = T)

#colnames(res_sum) = c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf",
#                      "gwas.exposure", "gwas.outcome", "PP.H4.threshold", "colocalized")

# save the causal snps conditional on H4 is true (PP.H4>=0.8)
fwrite(causal_snp, paste0("coloc/200kb/snp_h4_pp/", exposure, "_", outcomes, "_causal_snps_200kb_cc_v2.txt.gz"), col.names=T, row.names=F, quote=F, sep="\t", append = F)

# save the 95% credible sets
fwrite(cs_snp, paste0("coloc/200kb/95cs/", exposure, "_", outcomes, "_95_CredibleSets_200kb_cc_v2.txt"), col.names=T, row.names=F, quote=F, sep="\t", append = F)


#-----------------------------------------
#
# Coloc between exposure and outcome gwas 
# Susie
#
# primary analysis: p12 = 1e-5
#-----------------------------------------

library(data.table)
library(dplyr)
library(coloc)
library(tidyr)

#rm(list=ls())

#install.packages("LDlinkR")
library(LDlinkR)

prot_tss = fread("ukb_prot_tss_949_04032024.txt") 

prot_tss = prot_tss %>% 
  mutate(Gene_name = folder_name) %>%
  separate(Gene_name, "Gene", sep = "_" )

prot_tss = prot_tss %>% 
  filter(Gene == "HEXIM1"|
           Gene =="CD46"|
           Gene =="MEGF10")

prot_tss  = prot_tss %>% 
  mutate(CHR = as.numeric(CHR)) %>%
  mutate(gene_start = as.numeric(gene_start)) %>%
  mutate(gene_end = as.numeric(gene_end))

#=====================================================================================================================================================

# because the LD matrix website does not support the submission of paralleled jobs 
# I used a for loop in R instead of the job array


# continuous outcomes
outcomes = c("WMH", "MD", "FA")

# dichotomous outcomes
#outcomes = c("EPVS", "LS_TOAST_MRI", "CMB")

for(j in 1:length(outcomes)) {
  for(i in 1:length(prot_tss$Gene)) { 
    tryCatch({
      
      exposure = prot_tss$folder_name[i]
      outcome = outcomes[j]
      cis_chr = prot_tss$CHR[i]
      gene = prot_tss$Gene[i]
      
      lbound = prot_tss$gene_start[i] - 200000
      lbound = ifelse(lbound < 0, 0, lbound)
      
      ubound = prot_tss$gene_end[i] + 200000
      
      print(paste0(cis_chr, ":", lbound, "-", ubound))
      
      inflam_file = fread(paste0("coloc/exposure_cis_aligned_file/", exposure, "_cis_for_coloc.txt.gz"))
      outcome_file  = fread(paste0("GWAS_Sumstats_Outcomes/hg38_version/by_chromosome/", outcome, "_", cis_chr, "_GWAS_Sumstats_GRCh38_cleaned.txt.gz"))
      
      inflam_200kb = inflam_file %>% 
        filter(BP >= lbound & BP <= ubound) 
      
      common_snps = inflam_200kb %>% 
        select(VarID) %>% 
        inner_join(outcome_file) %>%
        select(VarID)
      
      snps = common_snps$VarID
      
      ld_mat = LDmatrix(snps, pop = "CEU", r2d = "r2", token = "23bb68a0b932",               
                        file = FALSE,
                        genome_build = "grch38")
      
      rownames(ld_mat) = ld_mat$RS_number
      ld_mat = as.matrix(ld_mat[ ,-1])
      
      snp_na_index = which(is.na(ld_mat[ ,1]))
      
      ld_mat_2 = ld_mat[-snp_na_index, -snp_na_index]
      #dim(ld_mat_2)
      
      common_snps_wLD = as.data.frame(colnames(ld_mat_2))
      colnames(common_snps_wLD) = "VarID"
      
      df_exp = inflam_200kb %>% 
        inner_join(common_snps_wLD) %>%
        mutate(VARBETA = SE^2)
      
      df_out = outcome_file %>% 
        inner_join(common_snps_wLD) %>%
        mutate(VARBETA = SE^2)
      
      df_exp_list=list(
        snp= df_exp$VarID, 
        beta= df_exp$BETA, 
        varbeta= df_exp$VARBETA, 
        position= df_exp$BP, 
        N = median(df_exp$N),
        type= "quant", 
        sdY= 1, 
        LD = ld_mat_2)
      
      df_out_list=list(
        snp= df_out$VarID, 
        beta= df_out$BETA, 
        varbeta= df_out$VARBETA, 
        position= df_out$BP, 
        N = median(df_out$N),
        type= "quant",              # if the outcome is lacunar stroke, then change this to "cc"
        sdY= 1, 
        LD = ld_mat_2)
      
      check_dataset(df_exp_list, req="LD")
      check_dataset(df_out_list, req="LD")
      
      S_exp=runsusie(df_exp_list)
      summary(S_exp)
      
      S_out=runsusie(df_out_list)
      summary(S_out)
      
      if(requireNamespace("susieR",quietly=TRUE)) {
        susie.res=coloc.susie(S_exp,S_out)
        print(susie.res$summary)
      }
      
      summ = susie.res$summary
      summ$exposure = gene
      summ$outcome = outcome
      
      #colnames(summ) = c("nsnps",     "hit1",      "hit2",      "PP.H0.abf", "PP.H1.abf", "PP.H2.abf",
      #                   "PP.H3.abf", "PP.H4.abf", "idx1",      "idx2",      "exposure",  "outcome"  )
      
      fwrite(summ, 
             paste0("coloc/200kb_susie/res_coloc_susie_200kb.txt"), col.names=F, row.names=F, quote=F, sep="\t", append = T)
      
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }}


