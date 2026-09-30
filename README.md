# EWAS of early PTB by sex in the PAGE study

## Project Overview

The Placental Abruption Genetic Epidemiology study, PAGE, is a case-control study of placental abruption and preterm birth based in Lima, Peru. For this project, we sought to characterize the saliva DNA methylation profile of early preterm newborns by conducting sex-stratified Epigenome-wide association analyses, performing differentially methylated region (DMR) analysis, and in silico functional exploration of the main signals obtained in the EWAS.
   
## Research Question
   
Do infants born early preterm (< 34 gestational weeks) differ in their epigenetic profile from those born at term, and can these differences be observed when stratifying by sex? What are the sex-specific biological pathways most impacted by early PTB?
   
## Data Source

The phenotypic dataset for this project is the `Phenotype_ePTB_PAGE_n284_share.csv` found in the `data` folder, which includes the exposure (`eptb_bin`), main covariates, batch variables (surrogate variables) and a complete list of predicted saliva cells (epithelial and immune cells) using EpiDISH (https://bioconductor.org/packages/devel/bioc/vignettes/EpiDISH/inst/doc/EpiDISH.html?utm_source=openai). Quality-controlled DNA methylation data obtained using the EPIC array is in the 'TO COMPLETE' file. 


## Base code for running methylation PCs (mPCs)

Methylation PCs were calculated using the 10K most variable normalized betas. The code for this analysis is: `mPCs_dataStructure_inspection_preEWAS.R` found in the `code` folder


## EWAS of ePTB by sex

Based on the distribution of mPCs, we identified that the main variability in the first PC was driven by infant sex. Thus, the main EWAS were stratified by this variable. The code for this analysis is the `TO-COMPLETE.R` found in the `code` folder. Results of the EWAS were further adjusted using bacon, and the code for this is 'TO-COMPLETE.R'. We visually summarized results using Manhattan plots, QQ-plots, and volcano plots.


## Analysis workflow
The data and code accompanying this project can be used in the following order, consistent with the order in which analyses were conducted:

   1. Calculating predicted cells from saliva DNAm using EpiDISH
   2. Calculating DNA methylation PCs
   3. EWAS of early PTB by sex (boys and girls separately)
   4. Bacon adjustment of EWAS by sex
   5. Differentially methylated regional analysis (DMRs) by sex using bacon-adjusted estimates
   6. In silico functional exploration of main EWAS results by sex:
         - Genomic enrichment analysis
         - Enrichment for biological pathways using gProfiler
         - Look-up of main signals in other EWAS
       

## Results

Results obtained from the sex-stratified EWAS of early PTB, summary plots, DMRs by sex, and results from the functional analysis can be found in the `results` folder.




