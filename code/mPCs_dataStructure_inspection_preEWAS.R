###########################################
# PAGE
# ePTB and Newborn DNAm
# Study: PAGE - R21 ePTB
# Task: Inspecting funderlying population structure using methylation PCs
# Author: Diana Juvinao-Quintero
# Date: 10-17-2025
###########################################

# Tutorial of usig ewaff package for EWAS
# https://htmlpreview.github.io/?https://github.com/perishky/ewaff/blob/master/docs/tutorial.html

# # open an interactive session in FASRC
# salloc -p sapphire --mem 50G -t 0-08:00
# 
# # Open SPACK environment
# cd /n/holylabs/LABS/gelaye_lab/Lab/software/spack
# . share/spack/setup-env.sh
# 
# # install R module and others
# module load cmake/3.30.3-fasrc01 # required for meffil installation
# module load R # open R (v 4.5.1)
# R

# load libraries
# All packages installed again in a local library (09/25/25). Selected core from Bristol (#63)
# install.packages(c("remotes", "ggplot2", "gridExtra", "dplyr", "gdata", "data.table", "psych", "car")) 
# install.packages("devtools")       # if needed
# devtools::install_github("perishky/meffil")
# install_github("perishky/ewaff")
library(remotes)
library(meffil)
library(remotes)
library(ewaff)
library(meffil)
library(ggplot2)
library(gridExtra)
library(data.table)
library(dplyr)
library(gdata)
library(limma)

#set wd
setwd("/n/holylabs/gelaye_lab/Users/djuvinaoquintero/DNAm_PAGE_newborns")

############################################
# load the data
############################################
load("norm.beta.pc20.Robj")
load("ePTB_R21_PAGE_1st_pull_n304_04.11.25.RData")
cells<-as.data.frame(fread("saliva_cellcounts_PAGE_ref-Middleton_etal_2021.txt"))
samplesheet<-read.csv('RQ-027059_Gelaye_EPIC_Sample_Sheet_Final_052925.csv', skip=7) # provided by the lab. Will use plate as a covariate

############################################
#prep betas and samplepheno
############################################
betas<-norm.beta
samplepheno<- page_r21_pulled

# merge cells and batches with samplepheno
cells <- cells %>% rename(immune=CD45pos,epithelial=large)
cells$Sample_Name<- as.character(cells$Sample_Name)
samplesheet$Sample_Name<-substr(samplesheet$Sample_Name, 1, 5)

samplepheno<-left_join(samplepheno, cells, by=c("studyid"="Sample_Name"))
samplepheno<-left_join(samplepheno, samplesheet, by=c("studyid"="Sample_Name"))

dim(samplepheno) # 304  61
dim(betas)       # 928585    284  

keep(betas, samplepheno, EWAS_name, sure=T)
identical(samplepheno$studyid, colnames(betas)) # F

# match sampleids across datasets
samplepheno<- samplepheno[match(colnames(betas), samplepheno$studyid), ]
identical(samplepheno$studyid, colnames(betas)) # T
dim(samplepheno) # 284
table(samplepheno$casestat, useNA = "always") # case=145, controls=139

############################################
#check pheno data and recode if needed 
############################################
##check EWAS phenotypes:
# ePTB (cases n=145, controls n=139)
samplepheno$eptb_bin<-recode(samplepheno$ePTB, "moderate PTB"="ePTB case") 
# exposure
str(samplepheno$eptb_bin)
table(samplepheno$eptb_bin)
levels(samplepheno$eptb_bin)

#covariates:
#maternal age
summary(samplepheno$A2age) 
str(samplepheno$A2age) #numeric
sum(is.na(samplepheno$A2age)) # 2 NAs

#gravidity (instead of parity as there are more NAs for this latter)
samplepheno$gravidity<-as.factor(samplepheno$gravidity)
str(samplepheno$gravidity)
table(samplepheno$gravidity, useNA = "always") # 3 NAs

#maternal smoking behaviour
samplepheno$smokpreg<-as.factor(samplepheno$smokpreg)
str(samplepheno$smokpreg)
table(samplepheno$smokpreg, useNA = "always")

#maternal education
samplepheno$education<-as.factor(samplepheno$education)
str(samplepheno$education) 
table(samplepheno$education, useNA = "always") # 1 NA

#maternal planned pregnancy
samplepheno$e2planned<-as.factor(samplepheno$e2planned)
str(samplepheno$e2planned) 
table(samplepheno$e2planned, useNA = "always") # 2 NAs

#maternal employment
samplepheno$employed<-as.factor(samplepheno$employed)
str(samplepheno$employed)
table(samplepheno$employed, useNA = "always") # 1 NA

#maternal pre-pegnancy BMI 
summary(samplepheno$prebmi) 
str(samplepheno$prebmi) #numeric
sum(is.na(samplepheno$prebmi)) # 5 NAs

# child sex
samplepheno$child_sex<-as.factor(samplepheno$child_sex)
str(samplepheno$child_sex)
table(samplepheno$child_sex, useNA = "always")  

#slide (Sentrix_ID)
## from normalization betas analysis, Slide was associated with methylation PC2 (derived from 10K most variable CpGs)
# batch.variable	principal.component	test	  p.value	  estimate	r2
# Slide		          PC2	              F-test	3.53e-86	42			  0.863
samplepheno$Sentrix_ID<-as.factor(samplepheno$Sentrix_ID)
str(samplepheno$Sentrix_ID)
table(samplepheno$Sentrix_ID, useNA = "always")

#cell-type estimates
# "immune", "epithelial"
summary(samplepheno %>% select(immune, epithelial)) 
str(samplepheno %>% select(immune, epithelial)) #numeric
sum(is.na(samplepheno %>% select(immune, epithelial))) # 0

##drop to complete cases
samplepheno<- samplepheno[complete.cases(samplepheno[,c("eptb_bin", "A2age", "gravidity", "smokpreg", "education", "employed","prebmi","child_sex", "e2planned", "immune", "epithelial")]), ]
dim(samplepheno) # 271 samples

table(samplepheno$eptb_bin) # 132 controls 139 cases

m <- match(samplepheno$studyid,colnames(betas))
betas <- betas[,m]

identical(samplepheno$studyid, colnames(betas)) # T

##############################################################
# exclude missing values in beta matrix and match datasets
##############################################################
cpgs_missing<- rownames(betas)[rowSums(is.na(betas)) > 0] # 824 CpGs
# these 824 cpgs are "nv-GRCh38" CpGs in the EPIC V2 array, which are uniquely defined and supported in GRCh38, not in GRCh37/hg19

# check the distribution of missing methylation values for these CpGs by columns
missing_betas<-betas[cpgs_missing, ]
sum_samples_cpgs_nas<-apply(missing_betas,2, function(x){sum(is.na(x))})
print(sum_samples_cpgs_nas) # all 271 samples have NAs for these 824 CpGs, exclude from analysis

# exclude these cpgs with NAs
betas<-betas[!rownames(betas) %in% cpgs_missing, ]
dim(betas) # 927761    271

# verify beta matrix has no NAs left
get_cpgs_with_missing <- function(beta_mat) {
  # Find CpGs (rows) with at least one missing value
  missing_cpgs_idx <- which(rowSums(is.na(beta_mat)) > 0)
  # If rownames are present, return them; otherwise, return row indices
  if (!is.null(rownames(beta_mat))) {
    cpg_ids <- rownames(beta_mat)[missing_cpgs_idx]
  } else {
    cpg_ids <- missing_cpgs_idx
  }
  return(cpg_ids)
}

get_cpgs_with_missing(betas) # 0

# check that samples still match between datasets
identical(samplepheno$studyid, colnames(betas)) # T
dim(betas) ; dim(samplepheno)

# clear wkspce
keep(samplepheno, betas, sure=T) 


####################################
# PCs from most variable probes
####################################
require(limma)
library(RColorBrewer)
library(matrixStats)
## select top 10K CpGs with largest variance
variable_betas_10K <- betas[order(-matrixStats::rowVars(betas))[1:10000],]
variable_betas_10K<-t(variable_betas_10K)
res.pca_10K <- prcomp(variable_betas_10K, scale = TRUE) #people are rows, CpGs are cols
summary(res.pca_10K)
pcs_10K<-res.pca_10K$x

## Generate screeplot
# Calculate proportion of variance explained for each PC
explained_var <- res.pca_10K$sdev^2 / sum(res.pca_10K$sdev^2)
cum_explained <- cumsum(explained_var)
tiff("PC_plots/screeplot_cumVar_methPCs_top10Kcgs.tiff",units="cm", width=12, height=8, res=300)
fviz_screeplot(res.pca_10K, addlabels=FALSE, ylim=c(0, 75), 
               choice="variance",
               xlab="Principal Component", 
               ylab="Cumulative Proportion of Variance",
               ncp=20,  # number of dimensions/pcs to plot
               geom="line", 
               geom.args = list(size = 0.5),  # adjust point size in the plot
               main="Variance Explained by top 20 mPCs",
               ggtheme = theme_classic())
dev.off()

pcs_10K<-as.data.frame(pcs_10K)
pcs_10K$studyid<-rownames(pcs_10K)
samplepheno<- left_join(samplepheno, pcs_10K[, c("studyid", "PC1", "PC2")], by="studyid")

####################################
# PC1 & PC2 by Technical Artifacts
####################################

# Sample Plate
tiff("PC_plots/norm_betas_10K_mostVar_probes_PC1&PC2_bySamplePlate.tiff",units="cm", width=13, height=8, res=300, pointsize=7)
plot(samplepheno$PC1, samplepheno$PC2, col=as.numeric(as.factor(samplepheno$Sample_Plate)),
     xlab="PC1", ylab="PC2", main="PC1 & PC2 ~ Sample Plate", pch=16)
legend("topright", legend=levels(as.factor(samplepheno$Sample_Plate)), col=1:length(levels(as.factor(samplepheno$Sample_Plate))), pch=16)
dev.off()

# Sample Well (multiple levels) 
tiff("PC_plots/norm_betas_10K_mostVar_probes_PC1&PC2_bySampleWell.tiff",units="cm", width=13, height=8, res=300, pointsize=7)
par(mar=c(5, 4, 4, 8)) 
plot(samplepheno$PC1, samplepheno$PC2, col=as.numeric(as.factor(samplepheno$Sample_Well)),
     xlab="PC1", ylab="PC2", main="PC1 & PC2 ~ Sample Well", pch=16)
# plot legend outside plot area
legend("topright", legend=levels(as.factor(samplepheno$Sample_Well)), col=1:length(levels(as.factor(samplepheno$Sample_Well))), pch=16, cex=0.4, inset = c(-0.1, 0), xpd=TRUE)
dev.off()

# Sentrix_ID (multiple levels)
tiff("PC_plots/norm_betas_10K_mostVar_probes_PC1&PC2_bySentrixID.tiff",units="cm", width=13, height=8, res=300, pointsize=7)
par(mar=c(5, 4, 4, 8)) 
plot(samplepheno$PC1, samplepheno$PC2, col=as.numeric(as.factor(samplepheno$Sentrix_ID)),
     xlab="PC1", ylab="PC2", main="PC1 & PC2 ~ Sentrix ID", pch=16)
# plot legend outside plot area
legend("topright", legend=levels(as.factor(samplepheno$Sentrix_ID)), col=1:length(levels(as.factor(samplepheno$Sentrix_ID))), pch=16, cex=0.4, inset = c(-0.1, 0), xpd=TRUE)
dev.off()

# Sentrix Position
tiff("PC_plots/norm_betas_10K_mostVar_probes_PC1&PC2_bySentrixPosition.tiff",units="cm", width=13, height=8, res=300, pointsize=7)
plot(samplepheno$PC1, samplepheno$PC2, col=as.numeric(as.factor(samplepheno$Sentrix_Position)),
     xlab="PC1", ylab="PC2", main="PC1 & PC2 ~ Sentrix Position", pch=16)
legend("topleft", legend=levels(as.factor(samplepheno$Sentrix_Position)), col=1:length(levels(as.factor(samplepheno$Sentrix_Position))), pch=16, cex=0.5)
dev.off()


####################################
# PC1 & PC2 by Biological variables
####################################
# brewer.pal.info # get palettes info
# pal <- brewer.pal(8,"Dark2") # choose a palette with as many colors as groups existent in grouping variable

# Case status (ePTB)
tiff("PC_plots/norm_betas_10K_mostVar_probes_PC1&PC2_byePTB.tiff",units="cm", width=13, height=8, res=300, pointsize=7)
plot(samplepheno$PC1, samplepheno$PC2, col=as.numeric(samplepheno$eptb_bin),
     xlab="PC1", ylab="PC2", main="PC1 & PC2 ~ Case status (ePTB)", pch=16)
legend("topright", legend=levels(as.factor(samplepheno$eptb_bin)), col=1:length(levels(as.factor(samplepheno$eptb_bin))), pch=16)
dev.off()


# Child sex 
tiff("PC_plots/norm_betas_10K_mostVar_probes_PC1&PC2_bychildSex.tiff",units="cm", width=13, height=8, res=300, pointsize=7)
plot(samplepheno$PC1, samplepheno$PC2, col=as.numeric(samplepheno$child_sex),
     xlab="PC1", ylab="PC2", main="PC1 & PC2 ~ child Sex", pch=16)
legend("topright", legend=levels(as.factor(samplepheno$child_sex)), col=1:length(levels(as.factor(samplepheno$child_sex))), pch=16)
dev.off()


# Gravity
tiff("PC_plots/norm_betas_10K_mostVar_probes_PC1&PC2_byGravidity.tiff",units="cm", width=13, height=8, res=300, pointsize=7)
plot(samplepheno$PC1, samplepheno$PC2, col=as.numeric(samplepheno$gravidity),
     xlab="PC1", ylab="PC2", main="PC1 & PC2 ~ Gravidity", pch=16)
legend("topright", legend=levels(as.factor(samplepheno$gravidity)), col=1:length(levels(as.factor(samplepheno$gravidity))), pch=16)
dev.off()


# Employment
tiff("PC_plots/norm_betas_10K_mostVar_probes_PC1&PC2_byEmployment.tiff",units="cm", width=13, height=8, res=300, pointsize=7)
plot(samplepheno$PC1, samplepheno$PC2, col=as.numeric(samplepheno$employed),
     xlab="PC1", ylab="PC2", main="PC1 & PC2 ~ Employment", pch=16)
legend("topright", legend=levels(as.factor(samplepheno$employed)), col=1:length(levels(as.factor(samplepheno$employed))), pch=16)
dev.off()


# Smoking
tiff("PC_plots/norm_betas_10K_mostVar_probes_PC1&PC2_bySmoking.tiff",units="cm", width=13, height=8, res=300, pointsize=7)
plot(samplepheno$PC1, samplepheno$PC2, col=as.numeric(samplepheno$smokpreg),
     xlab="PC1", ylab="PC2", main="PC1 & PC2 ~ Smoking", pch=16)
legend("topright", legend=levels(as.factor(samplepheno$smokpreg)), col=1:length(levels(as.factor(samplepheno$smokpreg))), pch=16)
dev.off()


# Planned pregnancy
tiff("PC_plots/norm_betas_10K_mostVar_probes_PC1&PC2_byPlannedPreg.tiff",units="cm", width=13, height=8, res=300, pointsize=7)
plot(samplepheno$PC1, samplepheno$PC2, col=as.numeric(samplepheno$e2planned),
     xlab="PC1", ylab="PC2", main="PC1 & PC2 ~ Planned Pregnancy", pch=16)
legend("topright", legend=levels(as.factor(samplepheno$e2planned)), col=1:length(levels(as.factor(samplepheno$e2planned))), pch=16)
dev.off()


## save dataset including PCs
save(samplepheno, file="ePTB_R21_PAGE_postProcessing_n271_2PCs_10.20.25.RData")

### END
