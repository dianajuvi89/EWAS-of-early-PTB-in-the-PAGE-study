# EWAS-ePTB-Aim2-project
Sharing code for EWAS-ePTB analyses Aim 2A


# Project Overview

The PAGE study is a case control study of placental abruption and preterm birth based in lima, Peru. For Aim 2A of this NICHD grant, we want to characterize the DNA methylation profile of early preterm infants, compare signals obtained with those of the continuous exposure, gestational age at delivery, and conduct functional analyses.
   
## Research Question
   
Does being born Early Preterm (< 32 gestational weeks) shapes differentially the methylome of newborn cases versus controls? What are the biological pathways most impacted by ePTB? are there sex-differences in these associations?
   
## Data Source

The main phenotype dataset for this project is the `Pheno_ePTB_PAGE_postBetasNormalizat_n284_CellsComplete.csv` found in the `data` folder, which includes the exposure (`eptb_bin`), main covariates, and a complete list of predicted saliva cells (epithelial and immune cells) using EpiDISH (https://bioconductor.org/packages/devel/bioc/vignettes/EpiDISH/inst/doc/EpiDISH.html?utm_source=openai).


# Base code for running methylation PCs (mPCs)

Methylation PCs were calculated using the 10K most variable normalized betas. 
The code for this analysis is: `mPCs_dataStructure_inspection_preEWAS.R` found in the `code` folder


# EWAS of ePTB by sex

Based on the distribution of mPCs, we identified that main variability in the first PC was driven by child sex. Thus, we decided to run the EWAS stratified by this variable. 

Code for this analysis is the `EWAS_ePTB_newborn_DNAm_cc-noEosinos_allCells_boys.R` and `EWAS_ePTB_newborn_DNAm_cc-allCells_girls.R` found in the `code` folder. The code includes plotting of main findings (Manhattan, QQ-plot, Volcano)

Note: the analysis for boys excludes adjustment for Eosinophils, as these immune cells were strongly correlated with the `eptb_bin` exposure in this dataset.


# Adjustment of EWASs using Bacon and Plots

This code conducts the adjustment of EWAS estimates (coef, se) using `bacon` and generates all plots based on bacon-corrected estimates. It also generates top sites derived from EWAS at p<`9.4e-08`

The code is the `Bacon-adj-topsites_EWAS_ePTB_cc-allCells_girls.R` and `Bacon-adj-topsites_EWAS_ePTB_cc-noEosinos_allCells_boys.R` and can be found in the `code` folder



# Preliminary results

All the results obtained in the main EWAS and after bacon adjustment will be stored in the `results` folder, including plots.




