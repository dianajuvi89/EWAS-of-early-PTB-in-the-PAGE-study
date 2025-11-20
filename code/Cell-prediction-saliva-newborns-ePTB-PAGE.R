#################################################################################
## Estimate cell-type fractions for surrogate saliva tissue using EpiDISH/HEpiDISH
## Running job in FAS-RC cluster in an interactive session
## Project: early PTB EWAS (R21 - Aim 2.A)
## Author: Diana Juvinao-Quintero
## Date: 11-18-2025
#################################################################################

# Example protocol to calculate cells here: https://bioconductor.org/packages/devel/bioc/vignettes/EpiDISH/inst/doc/EpiDISH.html?utm_source=openai

# Reference
# Zheng, Shijie C, Charles E Breeze, Stephan Beck, and Andrew E Teschendorff. 2018. “Identification of differentially methylated cell-types in Epigenome-Wide Association Studies.” Nature Methods 15 (12): 1059–66.

# Two steps of cell's estimation:
# The first step calculates the big cell fractions based on Zheng et al. 2018 (Epithelial, Fibroblasts, Immune cells) and the "RPC" method (Robust Partial Correlation). The second step calculates detail cell fractions for one of the 3 main cell types. Normally, fractions for immune cells (i.e., B, NK, CD8T, CD4T, Mono, Neutro, Eosino).  Epithelial and ICs tend to be the most prevalent cell types in saliva. 


# request memory
salloc -p sapphire --mem 40G -t 0-08:00 

# Open SPACK environment
cd /n/holylabs/LABS/gelaye_lab/Lab/software/spack
. share/spack/setup-env.sh

# install R module
module load R 
# open R (version 4.5.1)
R

#' remove all objects from your workspace
rm(list=ls())
gc()

# path to source file required to install the packages
source("/n/holylabs/LABS/gelaye_lab/Users/djuvinaoquintero/DNAm_PrOMIS-Kids/Data_QC/EPIC_QC/install_needed_packages.R") # PATH TO SOURCE FILE, this is an example path, change accordingly

#' Install packages if not already installed
need_pkgs <- c("FlowSorted.Blood.EPIC", "EpiDISH")
install_bioconductor_pkgs(pkgs = need_pkgs)

#' Run test
check_installed(pkgs = need_pkgs)

#' Load all packages, if needed
load_pkgs <- c("FlowSorted.Blood.EPIC", "EpiDISH", "data.table", "tibble", "feather", "dplyr")
lapply(load_pkgs, require, character.only = TRUE)


# set WD
setwd("/n/holylabs/LABS/gelaye_lab/Users/djuvinaoquintero/DNAm_PAGE_newborns")

#' This function will calculate the cell proportions and combine it with the phenotype information
#' input: beta values before Combat, phenotype file, reference cell proportions,
#' xcol- column name on which the data should be merged
#' ycol- column name
#' Output: Phenotye information with cell types
estimate_cell_proportion <- function(betavals, phenotyps, cell_prop_ref, xcol, ycol ){
  
  ## Calculate cell types using RPC method
  RPC <- epidish(betavals, as.matrix(cell_prop_ref), method = "RPC")
  
  cellTypes <- as.data.frame(RPC$estF) #RPC count estimates
  
  phen <- merge(phenotyps, cellTypes, by.x = xcol, by.y = ycol, all.x = T)
  
  return(phen)
}

## Load beta values (not ComBAT adjusted beta values, just normalized beta values)
#==============================================================================================
load("norm.beta.pc20.Robj")
beta <- norm.beta
rm(norm.beta)
dim(beta) # 928585  284

## Load ref EpiFibIC dataset (ref for Epithelial, Fibroblast and Immune Cell groups)
#==============================================================================================
data(centEpiFibIC.m)
ref<- centEpiFibIC.m
dim(ref) # 716  3

# check overlap between rownames in beta matrix and in ref cell matrix 
#==============================================================================================
table(rownames(beta) %in% rownames(centEpiFibIC.m)) # T 657 --> 657/716 DMCs for cells were found in our dataset

## Load the phenotype file
#==============================================================================================
load("ePTB_R21_PAGE_1st_pull_n304_04.11.25.RData")
pheno<- page_r21_pulled %>% mutate(eptb_bin=recode(ePTB, "moderate PTB"="ePTB case"))
# # save pheno as csv file
# write.csv(pheno, file="pheno_ePTB_R21_PAGE_1st_pull_n304_04.11.25.csv", row.names = FALSE)

## exclude samples that did not surpass QC
#==============================================================================================
library(dplyr)
# here, match two datasets as the normalized beta matrix has already excluded problematic samples based on the meffil protocol
pheno<- pheno %>% filter(studyid %in% colnames(beta))
pheno<- pheno[match(colnames(beta), pheno$studyid), ]

# make sure colnames in beta match with pheno
all(colnames(beta)==pheno$studyid) # T
identical(colnames(beta), pheno$studyid) # T

dim(pheno) # 284
dim(beta) # 284


## Run Cell prediction function
#==============================================================================================
# these are only fractions of Epi, IC, and Fibrinogen cells to see which cell type is most predominant
output <- estimate_cell_proportion(betavals = beta, phenotyps = pheno, cell_prop_ref = ref, xcol = "studyid", ycol = 'row.names')
# Save the phenotype file with cell types
# write.csv(output, file = "Pheno_ePTB_PAGE_postBetasNormalizat_n284_3CellFractions.csv", row.names = FALSE) 


# check the proportion of cell groups in the sample
#==============================================================================================
## For 95% of the samples, the distribution of cells was: 
quantile(output$Epi, 0.95, na.rm=T) # 95% of the samples had 0.9861381  or 98% of Epithelial cells
quantile(output$Fib, 0.95, na.rm=T) # 95% of the samples had 0.00       or 0.0% of Fibroblast cells
quantile(output$IC, 0.95, na.rm=T)  # 95% of the samples had 0.3775772  or 38% of Immune cells  --> highest cell group. Adjust for the subtypes of this cellular group  in the EWAS


# Estimate subtypes of ICs (second most prominent group of cells in saliva)
#==============================================================================================
#' This function will calculate the cell subtypes and combine it with the phenotype information
estimate_cell_proportion2 <- function(betavals, phenotyps, cell_prop_ref, cell_prop_ref2, xcol, ycol ){
  
  ## Calculate cell subtypes for ICs using RPC method in a two-step framework 
  ### h.CT.idx = 3 means that column 3 (i.e., ICs) from the fraction of cells estimated in first step, will be used in the second          independent reference data to calculate cell subtypes.
  RPC <- hepidish(betavals, as.matrix(cell_prop_ref), as.matrix(cell_prop_ref2), h.CT.idx = 3, method = "RPC", maxit=500) 
  # increase "maxit iterations". first round did not converse in 50 steps
  
  cellTypes <- as.data.frame(RPC) #RPC count estimates
  
  phen <- merge(phenotyps, cellTypes, by.x = xcol, by.y = ycol, all.x = T)
  
  return(phen)
}

# reference dataset for ICs in general tissues (i.e., not restricted to blood)
#==============================================================================================
data(centBloodSub.m)
ref2<-centBloodSub.m
## run the function
output2 <- estimate_cell_proportion2(betavals = beta, phenotyps = pheno, cell_prop_ref = ref, cell_prop_ref2=ref2,  xcol = "studyid", ycol = 'row.names' )

## save pheno with cells and subtypes
#==============================================================================================
write.csv(output2, file = "Pheno_ePTB_PAGE_postBetasNormalizat_n284_CellsComplete.csv", row.names = FALSE) # Save the phenotype file with cell types


## Open R in desktop and generate plot of cell proportions
#==============================================================================================
prelimCells<-read.csv("Pheno_ePTB_PAGE_postBetasNormalizat_n284_3CellFractions.csv")
cells<- read.csv("Pheno_ePTB_PAGE_postBetasNormalizat_n284_CellsComplete.csv")

# We can easily check the inferred fractions with boxplots. 
library(RColorBrewer)
pal <- brewer.pal(12,"Paired") 

## Fractions of the Three Major cell groups (Fib, Epi, ICs)
label=c("Epi", "Fib", "ICs")
tiff("Distribution_Epi-Fib-IC_cells_Zheng2018_saliva_ePTB_newborns_PAGE.tiff",units="cm", width=20, height=15, res=300, pointsize=7)
boxplot(prelimCells$Epi, prelimCells$Fib, prelimCells$IC, col=pal,
        ylab="Cellular Proportions", 
        main="Distribution Epithelial, Fibroblast and Immune Cells using EpiDISH \n Saliva samples of PAGE newborns (n=284)",
        names=label)
dev.off()

## Fractions of the Three Major cell groups (Fib, Epi, ICs) + IC subtypes
label=c("Epi", "Fib", "B", "NK", "CD4T", "CD8T", "Mono", "Neut", "Eosino")
tiff("Distribution_predCellsICsubtypes_saliva_ePTB_newborns_PAGE.tiff",units="cm", width=20, height=15, res=300, pointsize=7)
boxplot(cells$Epi, cells$Fib, cells$B, cells$NK, cells$CD4T, cells$CD8T, cells$Mono, cells$Neutro, cells$Eosino, col=pal,
        ylab="Cellular Proportions", 
        main="Distribution of Epithelial, Fibroblast and IC subtypes using EpiDISH \n Saliva samples of PAGE newborns (n=284)",
        names=label)
dev.off()


## Fraction of Immune Cell subtypes only
label=c("B", "NK", "CD4T", "CD8T", "Mono", "Neut", "Eosino")
tiff("Distribution_predICsSubtypes_saliva_ePTB_newborns_PAGE.tiff",units="cm", width=20, height=15, res=300, pointsize=7)
boxplot(cells$B, cells$NK, cells$CD4T, cells$CD8T, cells$Mono, cells$Neutro, cells$Eosino, col=pal,
        ylab="Cellular Proportions", 
        main="Distribution of Immune Cells subtypes using HEpiDISH \n Saliva samples of PAGE newborns (n=284)",
        names=label)
dev.off()


## END