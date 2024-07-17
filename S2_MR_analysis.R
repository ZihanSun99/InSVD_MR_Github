

#--------------------------------------
#
# Mendelian Randomization Analysis
#
# By Zihan Sun 
#
#--------------------------------------

#-----------------------------------------------------------------------------------------------------
# Chapter 1. Lift over outcome GWAS sumstats from hg 19 to hg 38 
# 
# LS -- GWAS Sumstats originally prepared by Alexy on hg38; no lift over needed 
# FA, MD, PSMD, WMH -- GWAS Sumstats originally prepared by Fatemeh on hg 19; lifted over to hg38 
# CMB -- GWAS Sumstats from Knol et al. Neurology paper on hg 19; merged with hg38 panel using rs id
# EPVS-WM -- GWAS Sumstats from Dupperron et al. Nature Medicine paper on hg 19; lifted over to hg 38
# 
#-----------------------------------------------------------------------------------------------------

library(data.table)
library(dplyr)
#rm(list=ls())

#BiocManager::install("liftOver")
library(liftOver)
library(rtracklayer)
library(GenomicRanges)


#' Liftover GWAS from hg19 to hg38
#'
#' @param GWAS GWAS with a minimum of 2 columns CHR and BP. If
#'   additional columns included these will be preserved.
#' @param path_to_chain Path to UCSC chain file for transformation from hg19 to
#'   hg38 coordinates.
#'
#' @return GWAS with hg38 coordinates
#' @export


liftover_hg19_to_hg38 <- function(GWAS, path_to_chain){
  
  # Import chain file
  hg19_to_hg38 <- rtracklayer::import.chain(path_to_chain)
  
  # If GWAS CHR column does not have "chr" in name, add to allow liftover
  if(!stringr::str_detect(GWAS$CHR[1], "chr")){
    
    GWAS <- GWAS %>%
      dplyr::mutate(CHR = stringr::str_c("chr", CHR))
    
  }
  
  # Convert GWAS to GRanges object
  GWAS_GR <-
    GenomicRanges::makeGRangesFromDataFrame(GWAS,
                                            keep.extra.columns = TRUE,
                                            ignore.strand = TRUE,
                                            seqinfo = NULL,
                                            seqnames.field = "CHR",
                                            start.field = "BP",
                                            end.field = "BP",
                                            starts.in.df.are.0based = FALSE)
  
  GWAS_hg38 <-
    rtracklayer::liftOver(GWAS_GR, hg19_to_hg38) %>%
    unlist() %>%
    as.data.frame() %>%
    dplyr::rename(CHR = seqnames,
                  BP = start) %>%
    dplyr::select(-end, -width, -strand) %>%
    dplyr::mutate(CHR = stringr::str_replace(CHR, "chr", ""))
  
  return(GWAS_hg38)
  
}

#=================================================================================================================================

#----------------------------
# FA, MD, PSMD, WMH - Fatemeh
#----------------------------

markers = c("FA", "MD", "PSMD", "WMH")

for(i in 1:length(markers)) {
  tryCatch({
    
    # load the file
    file = fread(paste("hpc-work/GWAS_Sumstats_Outcomes/hg19_version/UKB_",markers[i],"_GWAS_EUR_May2023.txt",sep=''))
    
    #dim(file) # 9309262      15
    #head(file)
    #range(file$A2FREQ) # A2FREQ is the EAF not the MAF
    
    # call the liftover function
    file.hg38 = liftover_hg19_to_hg38(file, "hpc-work/tools/hg19ToHg38.over.chain")
    dim(file.hg38)
    
    # post-processing / cleaning after the lift over
    
    file.hg38 <- file.hg38 %>%
      filter(!(CHR == "X")) %>%
      filter(!(CHR == "Y")) %>%
      filter(!is.na(CHR)) %>%
      filter(!is.na(BP)) %>%
      mutate(CHR = as.numeric(CHR)) %>%
      mutate(BP = as.numeric(BP))
    
    write.table(file.hg38, 
                gzfile(paste("hpc-work/GWAS_Sumstats_Outcomes/hg38_version/",markers[i],"_GWAS_Sumstats_GRCh38.txt.gz",sep='')),
                col.names=T, row.names=F, quote=F, sep="\t")
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}




# check the numbers before and after the conversion

#--------------------
# MARKER BEFORE AFTER
#--------------------
# FA  9309262 9302280
# MD  9309252 9302270 
# PSMD 9309266 9302284 
# WMH 9309293 9302311
#====================

#=================================================================================================================================

#--------------
# EPVS - WM 
#--------------

# wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90244001-GCST90245000/GCST90244151/GCST90244151_buildGRCh37.tsv.gz

# load the file
file = fread("GWAS_Sumstats_Outcomes/GCST90244151_buildGRCh37.tsv.gz")
dim(file) # 8783191      11
head(file)

# have to rename the column names to be compatible with the function!
file = as.data.frame(file)
file = file %>% 
  mutate(CHR = chromosome, BP = base_pair_location)

# call the liftover function
file.hg38 = liftover_hg19_to_hg38(file, "hpc-work/tools/hg19ToHg38.over.chain")

# post-processing / cleaning after the lift over

file.hg38 <- file.hg38 %>%
  mutate(CHR = stringr::str_replace(CHR, "chr", "")) %>%
  filter(!(CHR == "X")) %>%
  mutate(CHR = as.numeric(CHR))

dim(file.hg38)      #8774491      13
head(file.hg38)

file.hg38 = file.hg38[ ,-c(3, 4)]  # remove the extra chromosome and base_pair_position columns from the hg19 coordinates

#CHR       BP effect_allele other_allele effect_allele_frequency    beta
#1  10  9958055             a            g                  0.5497 -0.0027
#2  10 98240868             a            g                  0.5712  0.0058
#3  10 98240888             a            c                  0.7917 -0.0018
#4  10 98242110             t            c                  0.0130  0.0064
#5  10 98242707             t            c                  0.9865  0.0174
#6  10 98243485             t            g                  0.8873  0.0076

#standard_error Zscore p_value N_EVENT     N
#1         0.0037 -0.417 0.67700    7668 30950
#2         0.0037  2.036 0.04175    9324 38598
#3         0.0044 -0.880 0.37910    9324 38598
#4         0.0168  0.574 0.56620    7626 30916
#5         0.0165  1.056 0.29120    7123 28492
#6         0.0058  1.461 0.14400    9244 38031


write.table(file.hg38, 
            gzfile("GWAS_Sumstats_Outcomes/hg38_version/EPVS_GWAS_Sumstats_GRCh38.txt.gz"),
            col.names=T, row.names=F, quote=F, sep="\t")


which(is.na(file.hg38$BP))
which(is.na(file.hg38$CHR))
class(file.hg38$BP)
class(file.hg38$CHR)

#file.hg38$BP = as.numeric(file.hg38$BP)

#--------------------------------------
# LS ~ European ancestry (~ 6000 cases)
#--------------------------------------

# wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014122/LS_GCST90014122_buildGRCh37.tsv

# load the file
file = fread("GWAS_Sumstats_Outcomes/original_GWAS_sumstats/LS_GCST90014122_buildGRCh37.tsv")
dim(file)  # 6932926       9
head(file)

# have to rename the column names to be compatible with the function!
file = as.data.frame(file)
file = file %>% 
  mutate(CHR = chromosome, BP = base_pair_location)

# call the liftover function
file.hg38 = liftover_hg19_to_hg38(file, "hpc-work/tools/hg19ToHg38.over.chain")

# post-processing / cleaning after the lift over

file.hg38 <- file.hg38 %>%
  mutate(CHR = stringr::str_replace(CHR, "chr", "")) %>%
  filter(!(CHR == "X")) %>%
  mutate(CHR = as.numeric(CHR))

dim(file.hg38)      #6929700      11
head(file.hg38)

file.hg38 = file.hg38[ ,-c(4, 5)]  # remove the extra chromosome and base_pair_position columns from the hg19 coordinates

#CHR       BP variant_id effect_allele other_allele    beta standard_error
#1  10 98240868  rs7899632             G            A -0.0381         0.0209
#2  10 98240888 rs61875309             C            A -0.0033         0.0257
#3  10 98243485 rs12258651             G            T -0.0300         0.0330
#4  10 98243547 rs72828461             G            A  0.0706         0.0563
#5  10 98244028  rs1359508             C            T  0.0405         0.0217
#6  10 98244603  rs1048754             G            A  0.0033         0.0257

#p_value NStudies
#1 0.06871       13
#2 0.89630       12
#3 0.36250       13
#4 0.20950       11
#5 0.06211       12
#6 0.89770       12

write.table(file.hg38, 
            "GWAS_Sumstats_Outcomes/hg38_version/LS_TOAST_MRI_GWAS_Sumstats_GRCh38.txt",
            col.names=T, row.names=F, quote=F, sep="\t")

which(is.na(file.hg38$BP))
which(is.na(file.hg38$CHR))
class(file.hg38$BP)
class(file.hg38$CHR)

#file.hg38$BP = as.numeric(file.hg38$BP)

# DONE with LIFT-OVER from hg 19 to hg 38
# Cheers!  :)

#--------------------------------------------------
# 
# compile the bim file from Chr 1 to Chr 22
#
#--------------------------------------------------

df = fread("genetic_resources/LDSC/GRCh38/plink_files/1000G.EUR.hg38.1.bim")

for (i in 2:22) {
  
  df.append = fread(paste("genetic_resources/LDSC/GRCh38/plink_files/1000G.EUR.hg38.", i, ".bim", sep = ""))
  df = rbind(df, df.append)
  
}


write.table(df, "hpc-work/tools/1000G.EUR.hg38.chr1-22.bim", col.names=T, row.names=F, quote=F, sep="\t")
# all SNPs contained in this bim (.txt) file have an rs id. 


#---------------------------------------------------------------------------
#
# Chapter 2: Data Cleaning of SVD Imaging GWAS SumStats
#
# *** with hg 38 ref panel ***
#
#---------------------------------------------------------------------------

library(data.table)
library(ieugwasr)
library(genetics.binaRies)
library(dplyr)

#rm(list=ls())

source("R_script/manual_data_cleaning/GWAS_effect_strand.R")

################################################
## Read in reference panel

## autosomes

RP_1KGP_auto = as.data.frame(fread("hpc-work/tools/1000G.EUR.hg38.chr1-22.bim"))
colnames(RP_1KGP_auto) = c("CHR","VarID","Dist","BP","A1.ref","A0.ref")

# A0.ref is the other allele on the Ref panel 
# A1.ref is the effect allele on the Ref panel 

dim(RP_1KGP_auto)     #9991229       6
head(RP_1KGP_auto)

#################################################
## Part I: Pre-processing 
## 
## FA, MD, PSMD, WMH SumStats --> Done :) 
#################################################

## Read in and filter

markers = c("FA", "MD", "WMH")

for(i in 1:length(markers)) {
  tryCatch({
    
    summ = as.data.frame(fread(paste("GWAS_Sumstats_1000G_cleaned/hg38_wo_1000G_ref/", markers[i], "_GWAS_Sumstats_GRCh38.txt.gz", sep = "")))
    dim(summ)
    # 9309293      15
    
    #head(summ)
    #CHR    BP         SNP A1 A2    A2FREQ     INFO     N TEST        BETA
    #1   1 15777 rs568149713  G  A 0.9827900 0.354422 37322  ADD -0.03385360
    #2   1 16949 rs199745162  C  A 0.9811040 0.463087 37322  ADD -0.01374830
    #3   1 18849 rs533090414  G  C 0.0247382 0.514534 37322  ADD  0.01375840
    #4   1 52238   rs2691277  G  T 0.0206770 0.379503 37322  ADD  0.02751220
    #5   1 55164   rs3091274  A  C 0.0157249 0.384441 37322  ADD  0.01415820
    #6   1 57292 rs201418760  T  C 0.9773210 0.363339 37322  ADD -0.00285035
    
    #SE    CHISQ    LOG10P EXTRA         P
    #1 0.0431592 0.615266 0.3637010    NA 0.4328117
    #2 0.0359411 0.146324 0.1536180    NA 0.7020726
    #3 0.0297891 0.213315 0.1909920    NA 0.6441811
    #4 0.0382717 0.516769 0.3258530    NA 0.4722229
    #5 0.0433612 0.106614 0.1284080    NA 0.7440327
    #6 0.0372762 0.005847 0.0273119    NA 0.9390487
    
    
    ## Merge with reference panel and check alleles -- autosomes
    
    # merge with ref - by BP and CHR
    
    file = inner_join(summ, RP_1KGP_auto, by=c("CHR","BP"))
    
    #-----------------------
    # filter out the indels
    #-----------------------
    
    file$dir_1 = effect_strand("A0.ref","A1.ref","A1","A2", file)
    file$dir_2 = ifelse((file$A1.ref==file$A2 & file$A0.ref==file$A1),  1,
                        ifelse((file$A1.ref==file$A1 & file$A0.ref==file$A2), -1, NA))
    
    file$dir  = ifelse(!is.na(file$dir_1), file$dir_1, file$dir_2) 
    file = file[!is.na(file$dir),]
    file$BETA = as.numeric(file$BETA) * file$dir
    
    dim(file) #8082654      22
    
    file$EAF = ifelse(file$dir==1,as.numeric(file$A2FREQ),1-as.numeric(file$A2FREQ))
    file$P = as.numeric(file[,"P"])
    #file$P[which(file$P<1e-256)] = 1e-256
    file = rename(file, rsid = SNP)
    file$SNP = paste(file$CHR, file$BP, file$A0.ref, file$A1.ref, sep = ":")
    
    # get columns
    file =  file[,c("SNP","VarID","CHR","BP","A0.ref","A1.ref","BETA","SE","P","N","EAF")]
    names(file) = c("SNP","VarID","CHR","BP","A0","A1","BETA","SE","P","N","EAF")
    
    file = file[order(file$BP),]
    file = file[order(file$CHR),]
    
    file = file %>%
      distinct(SNP, .keep_all = T)
    
    head(file)
    dim(file) # 8082654      12
    
    write.table(file,
                gzfile(paste("GWAS_Sumstats_1000G_cleaned/", markers[i], "_GWAS_Sumstats_GRCh38_cleaned.txt.gz",sep='')),
                col.names=T, row.names=F, quote=F, sep="\t")
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


#################################################
## Part II: Pre-processing 
## 
## MRI-CONFIRMED! 
## LS SumStats --> Done :)
#################################################

## Read in and filter MRI-confirmed LS (Alexy's data)

markers = "LS"

for(i in 1:length(markers)) {
  tryCatch({
    
    summ = as.data.frame(fread(paste("GWAS_Sumstats_1000G_cleaned/hg38_wo_1000G_ref/", markers[i], "_GWAS_Sumstats_GRCh38.txt", sep = "")))
    dim(summ)
    # 4571739      13
    
    colnames(summ) = c("CHR", "BP", "ID", "A0", "A1", "EAF", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")
    
    
    #head(summ)
    
    #CHR      BP               ID A0 A1       EAF     N TEST        BETA        SE
    #1   1 1121358 chr1:1185978:A:C  A  C 0.0931093 23118  ADD -0.00619287 0.0583386
    #2   1 1121480 chr1:1186100:T:C  T  C 0.0947050 23135  ADD -0.00152192 0.0575450
    #3   1 1121657 chr1:1186277:T:C  T  C 0.1725020 23136  ADD  0.02547990 0.0421238
    #4   1 1121715 chr1:1186335:C:T  C  T 0.0941121 23132  ADD  0.00467494 0.0576510
    #5   1 1122283 chr1:1186903:T:G  T  G 0.0940648 23133  ADD  0.00735531 0.0577185
    #6   1 1122468 chr1:1187088:T:C  T  C 0.1680070 23133  ADD  0.02214700 0.0426347
    
    #CHISQ     LOG10P EXTRA
    #1 0.01126860 0.03836050    NA
    #2 0.00069947 0.00926149    NA
    #3 0.36588000 0.26339800    NA
    #4 0.00657565 0.02901650    NA
    #5 0.01623950 0.04643520    NA
    #6 0.26984000 0.21936600    NA
    
    
    ## Merge with reference panel and check alleles -- autosomes
    
    # merge with ref - by BP and CHR
    
    file = inner_join(summ, RP_1KGP_auto, by=c("CHR","BP"))
    
    #-----------------------
    # filter out the indels
    #-----------------------
    
    file$dir_1 = effect_strand("A0.ref","A1.ref","A0","A1", file)
    file$dir_2 = ifelse((file$A1.ref==file$A1 & file$A0.ref==file$A0),  1,
                        ifelse((file$A1.ref==file$A0 & file$A0.ref==file$A1), -1, NA))
    
    file$dir  = ifelse(!is.na(file$dir_1), file$dir_1, file$dir_2) 
    file = file[!is.na(file$dir),]
    file$BETA = as.numeric(file$BETA) * file$dir
    
    dim(file) #4184831      20
    
    file$EAF = ifelse(file$dir==1,as.numeric(file$EAF),1-as.numeric(file$EAF))
    
    file$P = 10^(-file$LOG10P)
    #file$P[which(file$P<1e-256)] = 1e-256
    
    file$SNP = paste(file$CHR, file$BP, file$A0.ref, file$A1.ref, sep = ":")
    
    # get columns
    file =  file[,c("SNP","VarID","CHR","BP","A0.ref","A1.ref","BETA","SE","P","N","EAF")]
    names(file) = c("SNP","VarID","CHR","BP","A0","A1","BETA","SE","P","N","EAF")
    
    file = file[order(file$BP),]
    file = file[order(file$CHR),]
    
    file = file %>%
      distinct(SNP, .keep_all = T)
    
    head(file)
    dim(file) # 4184828      11
    
    #file = file %>% 
    #  mutate(VarID = paste(CHR, BP, sep = ":"))
    
    write.table(file,
                gzfile(paste("GWAS_Sumstats_1000G_cleaned/", markers[i], "_GWAS_Sumstats_GRCh38_cleaned.txt.gz",sep='')),
                col.names=T, row.names=F, quote=F, sep="\t")
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


#################################################
## Part III: Pre-processing 
## 
## CMB SumStats --> Done :)
## 
#################################################


summ = as.data.frame(fread("datasets/CMB-GWAS-2020/CHARGE-UKB-CMB-E1-P1-METAANALYSIS-noGC-indels-20190312_1.TBL.postmeta"))
dim(summ)
#6808110      19
head(summ)

library(ieugwasr)
library(genetics.binaRies)
#rm(list=ls())

source("R_script/manual_data_cleaning/GWAS_effect_strand.R")

################################################
## Read in reference panel

## autosomes

RP_1KGP_auto = as.data.frame(fread("hpc-work/tools/1000G.EUR.hg38.chr1-22.bim"))
colnames(RP_1KGP_auto) = c("CHR","VarID","Dist","BP","A1.ref","A0.ref")

# A0.ref is the other allele on the Ref panel 
# A1.ref is the effect allele on the Ref panel 

dim(RP_1KGP_auto)     #9991229       6
head(RP_1KGP_auto)

# !!! according to the Knol et al paper, Allele 1 is the effect allele and Allele 2 is the other allele. !!!

## Merge with reference panel and check alleles -- autosomes

# merge with ref - by BP and CHR

summ$Allele1 = toupper(summ$Allele1)
summ$Allele2 = toupper(summ$Allele2)

file = merge(summ, RP_1KGP_auto, by.x="MarkerName", by.y="VarID") 

#-----------------------
# filter out the indels
#-----------------------

file$dir_1 = effect_strand("A0.ref","A1.ref","Allele2","Allele1", file)
file$dir_2 = ifelse((file$A1.ref==file$Allele1 & file$A0.ref==file$Allele2),  1,
                    ifelse((file$A1.ref==file$Allele2 & file$A0.ref==file$Allele1), -1, NA))

file$dir  = ifelse(!is.na(file$dir_1), file$dir_1, file$dir_2) 
file = file[!is.na(file$dir),]
file$BETA = as.numeric(file$Effect) * file$dir

dim(file) #6578446      28

file$EAF = ifelse(file$dir==1,as.numeric(file$Freq1),1-as.numeric(file$Freq1))
file$P = as.numeric(file[,"P_value"])
#file$P[which(file$P<1e-256)] = 1e-256
file = rename(file, rsid = MarkerName)

# get columns
file =  file[,c("rsid","CHR.y","BP.y","A0.ref","A1.ref","BETA","StdErr","P","N", "Case_N")]
names(file) = c("VarID","CHR","BP","A0","A1","BETA","SE","P","N","Case_N")

file$SNP = paste(file$CHR, file$BP, file$A0, file$A1, sep = ":")

#head(file)

file = file[order(file$BP),]
file = file[order(file$CHR),]

file = file %>%
  distinct(SNP, .keep_all = T)

head(file)
dim(file) # 6578446      11

write.table(file,
            gzfile(paste("CMB_GWAS_Sumstats_GRCh38_cleaned_v2.txt.gz",sep='')),
            col.names=T, row.names=F, quote=F, sep="\t")

#################################################
## Part IV: Pre-processing 
## 
## EPVS SumStats --> Done :)
#################################################

## Read in and filter

summ = as.data.frame(fread(paste("GWAS_Sumstats_1000G_cleaned/hg38_wo_1000G_ref/EPVS_GWAS_Sumstats_GRCh38.txt.gz", sep = "")))

dim(summ)
# 8774491      11

#head(summ)
#CHR       BP effect_allele other_allele effect_allele_frequency    beta
#1  10  9958055             a            g                  0.5497 -0.0027
#2  10 98240868             a            g                  0.5712  0.0058
#3  10 98240888             a            c                  0.7917 -0.0018
#4  10 98242110             t            c                  0.0130  0.0064
#5  10 98242707             t            c                  0.9865  0.0174
#6  10 98243485             t            g                  0.8873  0.0076

#standard_error Zscore p_value N_EVENT     N
#1         0.0037 -0.417 0.67700    7668 30950
#2         0.0037  2.036 0.04175    9324 38598
#3         0.0044 -0.880 0.37910    9324 38598
#4         0.0168  0.574 0.56620    7626 30916
#5         0.0165  1.056 0.29120    7123 28492
#6         0.0058  1.461 0.14400    9244 38031

#summ = rename(summ, CHR = chromosome, BP = base_pair_location)

summ$effect_allele = toupper(summ$effect_allele)
summ$other_allele = toupper(summ$other_allele)

## Merge with reference panel and check alleles -- autosomes

# merge with ref - by BP and CHR

file = inner_join(summ, RP_1KGP_auto, by=c("CHR","BP"))


#-----------------------
# filter out the indels
#-----------------------

file$dir_1 = effect_strand("A0.ref","A1.ref","other_allele","effect_allele", file)

file$dir_2 = ifelse((file$A1.ref==file$effect_allele & file$A0.ref==file$other_allele),  1,
                    ifelse((file$A1.ref==file$other_allele & file$A0.ref==file$effect_allele), -1, NA))

file$dir  = ifelse(!is.na(file$dir_1), file$dir_1, file$dir_2) 
file = file[!is.na(file$dir),]
file$beta = as.numeric(file$beta) * file$dir

dim(file) #8267247      18

#file$EAF = ifelse(file$dir==1,as.numeric(file$effect_allele_frequency),1-as.numeric(file$effect_allele_frequency))
file$P = as.numeric(file[,"p_value"])
#file$P[which(file$P<1e-256)] = 1e-256
#file = rename(file, rsid = SNP)
file$SNP = paste(file$CHR, file$BP, file$A0.ref, file$A1.ref, sep = ":")

# get columns
file =  file[,c("SNP","VarID","CHR","BP","A0.ref","A1.ref","beta","standard_error","P","N","N_EVENT")]
names(file) = c("SNP","VarID","CHR","BP","A0","A1","BETA","SE","P","N","N_EVENT")

file = file[order(file$BP),]
file = file[order(file$CHR),]

file = file %>%
  distinct(SNP, .keep_all = T)

head(file)
dim(file) # 8434364      11
            

#file = file %>% 
#  mutate(VarID = paste(CHR, BP, sep = ":"))

write.table(file,
            gzfile(paste("GWAS_Sumstats_1000G_cleaned/EPVS_GWAS_Sumstats_GRCh38_cleaned.txt.gz",sep='')),
            col.names=T, row.names=F, quote=F, sep="\t")

####################################################
## Part V: Pre-processing 
## 
## LS (TOAST + MRI-confirmed) SumStats --> Done :)
####################################################

## Read in and filter
summ = fread("GWAS_Sumstats_1000G_cleaned/hg38_wo_1000G_ref/LS_TOAST_MRI_GWAS_Sumstats_GRCh38.txt")
dim(summ)
# 6929700       9

#head(summ)
#CHR       BP variant_id effect_allele other_allele    beta standard_error
#1:  10 98240868  rs7899632             G            A -0.0381         0.0209
#2:  10 98240888 rs61875309             C            A -0.0033         0.0257
#3:  10 98243485 rs12258651             G            T -0.0300         0.0330
#4:  10 98243547 rs72828461             G            A  0.0706         0.0563
#5:  10 98244028  rs1359508             C            T  0.0405         0.0217
#6:  10 98244603  rs1048754             G            A  0.0033         0.0257
#p_value NStudies
#1: 0.06871       13
#2: 0.89630       12
#3: 0.36250       13
#4: 0.20950       11
#5: 0.06211       12
#6: 0.89770       12


## Merge with reference panel and check alleles -- autosomes

# merge with ref - by BP and CHR

file = inner_join(summ, RP_1KGP_auto, by=c("CHR","BP"))

#-----------------------
# filter out the indels
#-----------------------

file$dir_1 = effect_strand("A0.ref","A1.ref","other_allele","effect_allele", file)

file$dir_2 = ifelse((file$A1.ref==file$effect_allele & file$A0.ref==file$other_allele),  1,
                    ifelse((file$A1.ref==file$other_allele & file$A0.ref==file$effect_allele), -1, NA))

file$dir  = ifelse(!is.na(file$dir_1), file$dir_1, file$dir_2) 
file = file[!is.na(file$dir),]
file$beta = as.numeric(file$beta) * file$dir

dim(file) #6919013      16

file$P = as.numeric(file$p_value)
file$SNP = paste(file$CHR, file$BP, file$A0.ref, file$A1.ref, sep = ":")

# get columns
file =  file[,c("SNP","VarID","CHR","BP","A0.ref","A1.ref","beta","standard_error","P")]
names(file) = c("SNP","VarID","CHR","BP","A0","A1","BETA","SE","P")

file = file[order(file$BP),]
file = file[order(file$CHR),]

file = file %>%
  distinct(SNP, .keep_all = T)

head(file)
dim(file) # 6918519       9

write.table(file,
            "GWAS_Sumstats_1000G_cleaned/LS_TOAST_MRI_GWAS_Sumstats_GRCh38_cleaned.txt",
            col.names=T, row.names=F, quote=F, sep="\t")


#---------------------------------------------------------------------------
#
# Chapter 3: Two-Sample Mendelian Randomization Analyses 
#
# *** with hg 38 ref panel ***
#
#---------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)
library("genetics.binaRies")
library(ieugwasr)
library(TwoSampleMR)

#============================

writeLines("Start Job Num")

job_id <- commandArgs(trailingOnly=TRUE)[1]
print(job_id)
job_id <- as.numeric(job_id)

#job_id = 1

#======================================================================================================================================================
# make a trait list with outcome strings 

outcomes = c("WMH", "MD", "FA", "CMB", "LS_TOAST_MRI", "EPVS")
prot_tss = fread("tools/ukb_prot_tss_949.txt") 

multi_tss_list = data.frame(folder_name = c(
  "FUT3_FUT5_P21217_Q11128_OID21013_v1_Neurology",
  "AMY1A_AMY1B_AMY1C_P0DUB6_P0DTE7_P0DTE8_OID30707_v1_Inflammation_II", 
  "EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology", 
  "CKMT1A_CKMT1B_P12532_OID20721_v1_Inflammation"))

prot_tss = anti_join(prot_tss, multi_tss_list, by = "folder_name")  
# 4 markers with multi-tss sites were removed; their MR were performed separately to account for the multiple tss (due to the assay targeting multiple proteins)
# for these proteins, if their target proteins were also in the same genetic region, the cis-region was defined as +/- 1 Mb according to the "earliest" gene starting site and the last gene end site

prot_tss_wo_XY = prot_tss %>% 
  filter(CHR != "X") %>%
  filter(CHR != "Y")
# double check that 23 markers whose cis pQTL on X or Y chromosomes were excluded. 

prot_tss_wo_XY  = prot_tss_wo_XY %>% 
  slice(rep(1:nrow(prot_tss_wo_XY), each = 6)) %>%     # "each" argument repeat 1st item in the list 6 times first and then goes on to the second item
  mutate(outcome = rep(outcomes, 945)) %>%   # without "each", the entire list is repeated (in the order that the list is given) 945 times
  mutate(CHR = as.numeric(CHR)) %>%
  mutate(gene_start = as.numeric(gene_start)) %>%
  mutate(gene_end = as.numeric(gene_end))

#======================================================================================================================================================

writeLines("Read-in Outcome GWAS Sumstats")

outcome_str = prot_tss_wo_XY$outcome[job_id]

trait2_df = fread(paste0("MR_04042024/outcome_gwas/", outcome_str, "_GWAS_Sumstats_GRCh38_cleaned.txt.gz"))

trait2_df = trait2_df[ , c("SNP","VarID","CHR","BP","A0","A1","BETA","SE","P")]

exposure = prot_tss_wo_XY$folder_name[job_id]

cis_chr = prot_tss_wo_XY$CHR[job_id]
chr_file = prot_tss_wo_XY$file_name[job_id]

trait1_df = fread(paste("UKB_PPP_Sumstats/", exposure, "/discovery_chr", cis_chr, "_", chr_file, ".gz", sep = ""))

#----------------------------
# process exposure sumstats
#----------------------------

lbound = prot_tss_wo_XY$gene_start[job_id] - 1000000
lbound = ifelse(lbound < 0, 0, lbound)

ubound = prot_tss_wo_XY$gene_end[job_id] + 1000000

trait1_cis_df = trait1_df %>%
  mutate(CHR = as.numeric(CHROM)) %>% 
  filter(GENPOS >= lbound & GENPOS <= ubound) %>%
  mutate(P = 10^(-LOG10P)) %>%
  rename(BP = GENPOS, 
         A1 = ALLELE1, 
         A0 = ALLELE0,
         EAF = A1FREQ) %>% 
  select(CHR, BP, A0, A1, BETA, P, SE, N, EAF)

#-----------------
# filter out MHC
#-----------------

mhc = trait1_cis_df %>% 
  filter(CHR == "6") %>% 
  filter(BP >= 25500000 & BP <= 34000000) 

trait1_cis_df = anti_join(trait1_cis_df, mhc)

#-------------------------
# import outcome sumstats
#-------------------------

merge_file = trait1_cis_df %>% 
  inner_join(trait2_df, by = c("CHR", "BP")) %>%
  distinct(CHR, BP, .keep_all = T)

comm_snp_exp_file = merge_file %>%
  select(CHR, BP, SNP, VarID, A1.x, A0.x, BETA.x, P.x, SE.x, N, EAF) %>%
  arrange(CHR, BP)

comm_snp_exp_file <- format_data(
  comm_snp_exp_file,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  snp_col = "VarID",
  beta_col = "BETA.x",
  se_col = "SE.x",
  effect_allele_col = "A1.x",
  other_allele_col = "A0.x",
  eaf_col = "EAF",
  pval_col = "P.x",
  samplesize_col = "N") %>%
  filter(!is.na(effect_allele.exposure) & !is.na(other_allele.exposure)) 

comm_snp_exp_file$id.exposure = exposure

comm_snp_out_file = merge_file %>%
  select(CHR, BP, SNP, VarID, A1.y, A0.y, BETA.y, SE.y, P.y)

comm_snp_out_file <- format_data(
  comm_snp_out_file,
  type = "outcome",
  snps = NULL,
  header = TRUE,
  snp_col = "VarID",
  beta_col = "BETA.y",
  se_col = "SE.y",
  effect_allele_col = "A1.y",
  other_allele_col = "A0.y",
  pval_col = "P.y")

comm_snp_out_file$id.outcome = outcome_str 

comm_snp_clump_use = comm_snp_exp_file %>%
  #select(SNP, pval.exposure) %>%   #here, SNP is the rsid set by the exposure_format function
  distinct(SNP, .keep_all = T) %>%
  rename(ID = SNP, P = pval.exposure)

fwrite(comm_snp_clump_use, paste0("MR_04042024/exposure_cis_file/", exposure, "_", outcome_str, "_clump_use.txt.gz"), col.names=T, row.names=F, quote=F, sep="\t", append = F)

#==================================================================================================================================================================================================================

writeLines("Clumping")
cmd = paste0("plink2 --pfile genetic_resources/LDSC/GRCh38/plink_files/pgen/1000G.EUR.hg38.", cis_chr, " --clump MR_04042024/exposure_cis_file/", exposure, "_", outcome_str, "_clump_use.txt.gz --clump-p1 5e-8 --clump-r2 0.01 --clump-kb 1000 --out MR_04042024/instruments/", exposure, "_", outcome_str, "_5e-8_instruments") 
print(cmd)

setwd("software/plink2")
system(cmd)

#==================================================================================================================================================================================================================

writeLines("Read clumped results")
clump_res = fread(paste0("MR_04042024/instruments/", exposure, "_", outcome_str, "_5e-8_instruments.clumps"))

instrument_snp_exp = clump_res %>%
  select(ID) %>%
  inner_join(comm_snp_clump_use, by = "ID") %>%
  rename(SNP = ID, p.exposure = P)


# get harmonized instrument snp data (with the outcome data included)
dat <- harmonise_data(instrument_snp_exp, comm_snp_out_file, action=3)
write.table(dat, paste0("MR_04042024/instruments/snp_", outcome_str, "_5e-8_04042024.txt"), col.names=F, row.names=F, quote=F, sep="\t", append = T)

#colnames(dat) = c("SNP",                  "effect_allele.exposure",  "other_allele.exposure", 
#                  "effect_allele.outcome",  "other_allele.outcome",    "beta.exposure",         
#                  "beta.outcome",           "eaf.exposure",           "eaf.outcome",           
#                  "remove",                 "palindromic",            "ambiguous",             
#                  "id.outcome",             "se.outcome",             "pval.outcome",          
#                  "samplesize.outcome",     "outcome",                "mr_keep.outcome",       
#                  "pval_origin.outcome",    "pval.exposure",          "se.exposure",           
#                  "samplesize.exposure",    "exposure",               "mr_keep.exposure",      
#                  "pval_origin.exposure",   "id.exposure",            "action",                
#                  "mr_keep")


writeLines("Get MR results")

# get the MR results 
res = mr(dat, method_list = c("mr_ivw_fe", "mr_egger_regression", "mr_wald_ratio"))
write.table(res, paste0("MR_04042024/raw_results_5e-8/res_raw_", outcome_str, "_5e-8_04042024.txt"), col.names=F, row.names=F, quote=F, sep="\t", append = T)
#colnames(res) = c("id.exposure","id.outcome","outcome","exposure","method","nsnp","b","se","pval")


# get heterogeneity
res_hetero= mr_heterogeneity(dat, method_list = c("mr_ivw_fe", "mr_egger_regression"))
write.table(res_hetero, paste0("MR_04042024/raw_results_5e-8/res_hetero_", outcome_str, "_04042024.txt"), col.names=F, row.names=F, quote=F, sep="\t", append = T)
#colnames(res_hetero) = c("id.exposure", "id.outcome", "outcome", "exposure", "method", "Q", "Q_df", "Q_pval")


# get pleiotropy test; i.e. egger intercept 
res_pleio = mr_pleiotropy_test(dat)
write.table(res_pleio, paste0("MR_04042024/raw_results_5e-8/res_pleio_", outcome_str, "_04042024.txt"), col.names=F, row.names=F, quote=F, sep="\t", append = T)
#colnames(res_pleio) = c("id.exposure", "id.outcome", "outcome", "exposure", "egger_intercept", "se", "pval")

# Done with the UKB proteomic data 
# Iceland 36K SomaScan data were analysed with the same methods. 
#============================================================================================================

