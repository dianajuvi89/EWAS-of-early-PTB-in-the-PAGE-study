--------------------------------------------------------------------------------
1-cell-prediction-saliva-newborns-ePTB-PAGE.R
--------------------------------------------------------------------------------
PURPOSE
    Calculate predicted cell proportions from saliva DNAm using EpiDISH

INPUTS
    - Normalised DNAm data
    - Phenotype data

OUTPUTS
    - Phenotype data with estimated cell proportions

DEPENDENCIES
    - FlowSorted.Blodd.EPIC r-package
    - EpiDISH r-package

--------------------------------------------------------------------------------
2-mPCs_dataStructure_inspection_preEWAS.R
--------------------------------------------------------------------------------
PURPOSE
    Calculate DNAm principal components (PCs)

INPUTS
    - Normalised DNAm data
    - Phenotype data with estimated cell proportions
    - EPIC samplesheet

OUTPUTS
    - PC diagnostic plots
    - phenotype data with cell proportions and PCs

--------------------------------------------------------------------------------
3-EWAS-sex-stratified.R
--------------------------------------------------------------------------------
PURPOSE
    Carry out epigenome wide association analysis adjusted for confounders and 
    structural artifacts (using surrogate variable analysis)

INPUTS
    - Normalised DNAm data
    - Phenotype data with estimated cell proportions
    - the name of the sex for which to carry out the analysis

OUTPUTS
    - EWAS results
    - Standard plots of EWAS results

--------------------------------------------------------------------------------
4-Bacon-adjustment.R
--------------------------------------------------------------------------------
PURPOSE
    Carry out a Bacon adjustment on the obtained EWAS results to reduce inflation

INPUTS
    - EWAS results for a single sex and the name of that sex
    - Normalised DNAm data
    - Phenotype data with estimated cell proportions

OUTPUTS
    - Bacon-adjusted EWAS results for a single sex
    - Figure containing plots of results

--------------------------------------------------------------------------------
5-DMR-analysis.R
--------------------------------------------------------------------------------
PURPOSE
    Identify differentially methylated regions (DMRs) using Bacon-adjusted EWAS
    results

INPUTS
    - Bacon-adjusted EWAS results for a single sex and the name of that sex
    - Normalised DNAm data

OUTPUTS
    - DMR analysis results
    - Figure of DMR analysis results

--------------------------------------------------------------------------------
6.1-Lussier-analysis.R
--------------------------------------------------------------------------------
PURPOSE
    Identify which types of genetic sites are overrepresented in the bacon-adjusted 
    EWAS results

INPUTS
    - Bacon-adjusted EWAS results for a single sex and the name of that sex
    - EPIC2 manifest containing DNAm array annotations
    - p-value and association strength coefficient thresholds used for inclusion

OUTPUTS
    - enrichment results for enhancers, DNA hypersensitive sites, and CpG islands
    - plot of these results on one figure

--------------------------------------------------------------------------------
6.2-enrichment-analysis.R
--------------------------------------------------------------------------------
PURPOSE

INPUTS

OUTPUTS