######################
3.1-Stratified-EWAS.R
######################
- Code that produces results for adjusted and unadjusted EWAS using a
restricted linear model with the plots for stratified data
- For boys, the adjusted model exclude Eosino cells due to perfect prediction
(all boys with non-zero Eosino cell count were a case)
- The programme takes 10 arguments, including 4 input and 6 output arguments

- Can be run using the command "snakemake -c1 EWAS-[SEX]"

Input:
- Robj of normalised beta values from processing methylation data
- The sample phenotypes with the cell proportions added
- samplesheet
- sex by which to stratify (Male or Female)
- whether this is the sensitivity analysis (removing individuals with Eosino cell counts above 0)

Output:
- PATH/TO/ewas_results.csv
- PATH/TO/ewas_unadj.html
- PATH/TO/ewas_adj.html
- PATH/TO/ewas_results.RDATA
- PATH/TO/ewas_plots_unadj.pdf
- PATH/TO/ewas_plots_adj.pdf

##########################
3.2-Bacon-Adjusted-EWAS.R
##########################
- Code that produces Bacon adjusted EWAS results for EWAS results from
previous programme stratified by sex
- The programme takes 10 arguments, icnluding 4 inputs and 6 outputs

- Can be run using the command "snakemake -c1 EWAS-BACON-[SEX]"

Input:
- EWAS results .csv file
- sex by which data was stratified (Male or Female)
- Robj of normalised beta values from processing methylation data
- The sample phenotypes with the cell proportions added

Output:
- PATH/TO/bacon_plots_unadj.pdf
- PATH/TO/bacon_plots_adj.pdf
- PATH/TO/DIR/bacon_results
- PATH/TO/bacon_allplots_unadj.pdf
- PATH/TO/bacon_allplots_adj.pdf
- PATH/TO/ewas_topsites.csv