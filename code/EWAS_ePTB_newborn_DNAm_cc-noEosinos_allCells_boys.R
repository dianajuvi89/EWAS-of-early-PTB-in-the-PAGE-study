###########################################
# PAGE
# ePTB and Newborn DNAm
# Study: PAGE - R21 ePTB
# Analysis only in boys (refining adjusted model and including detailed cells)
# Author: Diana Juvinao-Quintero
# Date: 11-18-2025
###########################################

# Tutorial of usig ewaff package for EWAS
# https://htmlpreview.github.io/?https://github.com/perishky/ewaff/blob/master/docs/tutorial.html

# # open an interactive session in FASRC
# salloc -p sapphire --mem 50G -t 0-08:00
# 
# Open SPACK environment
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

#set wd
setwd("/n/holylabs/gelaye_lab/Users/djuvinaoquintero/DNAm_PAGE_newborns")

############################################
# load the data
############################################
load("norm.beta.pc20.Robj")
samplepheno<-read.csv("Pheno_ePTB_PAGE_postBetasNormalizat_n284_CellsComplete.csv")
samplesheet<-read.csv('RQ-027059_Gelaye_EPIC_Sample_Sheet_Final_052925.csv', skip=7) # provided by the lab. Will use plate as a covariate

#set EWAS name
EWAS_name<-"ePTB_boys_noNAbetas_noEosinos_allCells_rlm"

############################################
#prep betas and samplepheno
############################################
betas<-norm.beta
rm(norm.beta)
# merge batches in samplesheet with samplepheno
samplesheet$Sample_Name<-substr(samplesheet$Sample_Name, 1, 5)
samplepheno$studyid<-as.character(samplepheno$studyid)
samplepheno<-left_join(samplepheno, samplesheet, by=c("studyid"="Sample_Name"))
dim(samplepheno) # 284  69
dim(betas)       # 928585    284      
# free wkspace
keep(betas, samplepheno, EWAS_name, sure=T)

## subset to Male's IDs
samplepheno<- samplepheno %>% filter(child_sex=="Male")
dim(samplepheno) # should be 149 Correct.
# match betas with samplepheno 
betas<- betas[, colnames(betas) %in% samplepheno$studyid]
dim(betas); dim(samplepheno)

# verify if samples match across datasets
identical(samplepheno$studyid, colnames(betas)) # T
all(samplepheno$studyid== colnames(betas)) # T
# how many cases and controls among girls only
table(samplepheno$casestat, useNA = "always") # ePTB-cases=59 controls=90
table(samplepheno$eptb_bin, useNA = "always") # ePTB-cases=59 controls=90

# # remove methylation outliers based on IQR method (skip)
# betas<- ewaff.handle.outliers(betas, method="iqr")[[1]] #sets methn elements that are outside 3*IQR to NA (per CpG)
# 
# # identified outliers
# which(is.na(betas), arr.ind=T)

############################################
#check pheno data and recode if needed 
############################################
##check EWAS phenotypes:

# exposure
str(samplepheno$eptb_bin)
table(samplepheno$eptb_bin)
levels(samplepheno$eptb_bin)

#covariates:
#maternal age
summary(samplepheno$A2age) 
str(samplepheno$A2age) #numeric
samplepheno %>% group_by(eptb_bin) %>% summarise(mean=mean(A2age, na.rm=T), sd=sd(A2age, na.rm=T))
sum(is.na(samplepheno$A2age)) # 1 NAs

#gravidity (instead of parity as there are more NAs for this latter)
samplepheno$gravidity<-as.factor(samplepheno$gravidity)
str(samplepheno$gravidity)
samplepheno %>% select(eptb_bin, gravidity) %>% table()
table(samplepheno$gravidity, useNA = "always") # 3 NAs

# #maternal smoking behaviour (not too frequent in or sample - Excluded)
# samplepheno$smokpreg<-as.factor(samplepheno$smokpreg)
# str(samplepheno$smokpreg)
# table(samplepheno$smokpreg, useNA = "always")

# #maternal education (Excluded - employment instead)
# samplepheno$education<-as.factor(samplepheno$education)
# str(samplepheno$education) 
# table(samplepheno$education, useNA = "always") #  NA
# 
# #maternal planned pregnancy (Excluded)
# samplepheno$e2planned<-as.factor(samplepheno$e2planned)
# str(samplepheno$e2planned) 
# table(samplepheno$e2planned, useNA = "always") #  NAs

#maternal employment
samplepheno$employed<-as.factor(samplepheno$employed)
str(samplepheno$employed)
samplepheno %>% select(eptb_bin, employed) %>% table() # controls have equal numbers of moms employed/unemployed
table(samplepheno$employed, useNA = "always") # 1 NA
#           employed
# eptb_bin  No Yes
# control   29  29
# ePTB case 41  49


#maternal pre-pegnancy BMI 
summary(samplepheno$prebmi) 
str(samplepheno$prebmi) #numeric
samplepheno %>% group_by(eptb_bin) %>% summarise(mean=mean(prebmi, na.rm=T), sd=sd(prebmi, na.rm=T)) # similar values between groups
sum(is.na(samplepheno$prebmi)) # 5 NAs
# eptb_bin     mean    sd
# <chr>       <dbl> <dbl>
# 1 control    25.6  5.03
# 2 ePTB case  25.2  5.02

# child sex
samplepheno$child_sex<-as.factor(samplepheno$child_sex)
str(samplepheno$child_sex)
table(samplepheno$child_sex, useNA = "always")

#slide (Sentrix_ID)
## from normalization betas analysis, Slide was associated with methylation PC2 (derived from 10K most variable CpGs)
# batch.variable	principal.component	test	  p.value	  estimate	r2
# Slide		          PC2	              F-test	3.53e-86	42			  0.863
# samplepheno$Sentrix_ID<-as.factor(samplepheno$Sentrix_ID)
# str(samplepheno$Sentrix_ID)
# table(samplepheno$Sentrix_ID, useNA = "always")

#cell-type estimates 
# "immune detailed cells" and "epithelial" (exclude Fibroblasts)
summary(samplepheno %>% select(Epi, B:Eosino)) 
str(samplepheno %>% select(Epi, B:Eosino)) #numeric
apply(samplepheno %>% filter(eptb_bin=="control") %>% select(Epi, B:Eosino), 2, function(x){summary(x)})
apply(samplepheno %>% filter(eptb_bin=="ePTB case") %>% select(Epi, B:Eosino), 2, function(x){summary(x)})
sum(is.na(samplepheno %>% select(Epi, B:Eosino))) # 0


# test correlation between covariates 
library(psych) # for correlation analysis
library(car) # calculate VIF in pairwise analysis

covars<- samplepheno %>% select(eptb_bin, A2age, gravidity, employed, prebmi, Epi, B:Eosino)
levels<-lapply(covars[, c("eptb_bin", "gravidity", "employed")], function(x) {print(levels(x))})
covars<-covars %>% mutate(eptb_bin=dplyr::recode(eptb_bin, "control"=0, "ePTB case"=1),
                          gravidity=dplyr::recode(gravidity, "primiparous"=0, "multiparous"=1),
                          # smokpreg=dplyr::recode(smokpreg, "No"=0, "Yes"=1),
                          # education=dplyr::recode(education, "> high school"=0, "<= high school"=1),
                          employed=dplyr::recode(employed, "Yes"=0, "No"=1)
                          # child_sex=dplyr::recode(child_sex, "Female"=0, "Male"=1),
                          # e2planned=dplyr::recode(e2planned, "Yes"=0, "No"=1)
                          )
# Moderately high corr between: mom-age~gravid=0.4 (p<0.001) 
cor_matrix <- round(cor(covars, use = "pairwise.complete.obs"), 3) 
cor_test <- corr.test(as.matrix(covars))$p
print(cor_matrix); print(cor_test)

# Assess multicolinearity using variance inflation factor (VIF > 5 or 10 = high multicolinarity)
require(car)
vif_mod <- glm(eptb_bin ~ A2age + gravidity + employed + prebmi + Epi + B + CD4T + CD8T + Mono +NK + Neutro, data = covars) # leave Eosinos aside as they will cause Errors --> multicollinear with ePTB status in boys analysis only
# error while running vif() -> only ePTB cases have 0.0 values for Eosinos and none of the controls (for girls is equal)
vif(vif_mod) # high VIF between ePTB and cells

##drop to complete cases
samplepheno<- samplepheno[complete.cases(samplepheno[,c("eptb_bin", "A2age", "gravidity", "employed","prebmi", "Epi", "B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino")]), ]
dim(samplepheno) # 139 69

table(samplepheno$eptb_bin) # 55 controls 84 cases

# match again betas and pheno
m <- match(samplepheno$studyid,colnames(betas))
betas <- betas[,m]
identical(samplepheno$studyid, colnames(betas)) # T

dim(betas)       # 928585    139
dim(samplepheno) # 139  69

##############################################################
# exclude missing values in beta matrix and match datasets
##############################################################
cpgs_missing<- rownames(betas)[rowSums(is.na(betas)) > 0] 
length(cpgs_missing) # 824 CpGs
# these 824 cpgs are "nv-GRCh38" CpGs in the EPIC V2 array, which are uniquely defined and supported in GRCh38, not in GRCh37/hg19

# check the distribution of missing methylation values for these CpGs by columns
missing_betas<-betas[cpgs_missing, ]
sum_samples_cpgs_nas<-apply(missing_betas,2, function(x){sum(is.na(x))})
length(sum_samples_cpgs_nas) # all 139 samples have NAs for these 824 CpGs, exclude from analysis

# exclude these cpgs with NAs
betas<-betas[!rownames(betas) %in% cpgs_missing, ]
dim(betas) # 927761   XX

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
all(samplepheno$studyid== colnames(betas)) # T
dim(betas) ;     # 927761    139
dim(samplepheno) # 139  69


# test beta matrix
# betas<- betas[1:1000,]

##############################################################
# evaluate rank issues in model.matrix before EWAS with RLM
##############################################################
# Step 1. Build the model matrix
# This expands factors into dummy variables and builds the exact matrix EWAFF will use.
mm <- model.matrix(~ eptb_bin + A2age + gravidity + employed + prebmi + Epi + B + NK + CD4T + CD8T + Mono + Neutro + Eosino, data = samplepheno)
# Step 2. Check rank
# If rank < number of columns → collinearity problem (will cause RLM to fail)
qr(mm)$rank # 13 independent columns
ncol(mm) # 14
# Step 3. Identify which predictors cause the problem
# Use alias() to see dependencies:
alias(glm(rnorm(nrow(mm)) ~ mm))
# Model :
# rnorm(nrow(mm)) ~ mm
# Complete :
#            (Intercept) mmeptb_binePTBcase mmA2age mmgravidityprimiparous
# mm(Intercept)  1           0                   0       0                    
# mmEosino       1           0                   0       0                    
#               mmemployedYes mmprebmi mmEpi mmB  mmNK mmCD4T mmCD8T mmMono
# mm(Intercept)  0             0        0     0   0    0      0      0    
# mmEosino       0             0       -1    -1  -1   -1     -1     -1    
# mmNeutro
# mm(Intercept)  0      
# mmEosino      -1

# Note: Eosino may be causing errors in the only-boys analysis. Try excluding this variable and re-test rank issues in model.matrix
mm2<-model.matrix(~ eptb_bin + A2age + gravidity + employed + prebmi + Epi + B + NK + CD4T + CD8T + Mono + Neutro, data = samplepheno)
qr(mm2)$rank # 13  correct!
ncol(mm2)    # 13 correct!
alias(glm(rnorm(nrow(mm2)) ~ mm2))

# clear wkspce
keep(samplepheno, betas, EWAS_name, sure=T) 

############################################
#Run EWAS using ewaff
############################################

ewas.ret <- ewaff.sites(methylation ~ eptb_bin,
                        variable.of.interest="eptb_bin",
                        methylation=betas,
                        #generate.confounders="sva",
                        #n.confounders =10,
                        data=samplepheno,
                        method="rlm",
                        family="gaussian")

ewas.ret.adj <- ewaff.sites(methylation ~ eptb_bin + A2age + gravidity + employed + prebmi + 
                              Epi + B + NK + CD4T + CD8T + Mono + Neutro, # 
                            variable.of.interest="eptb_bin",
                            methylation=betas,
                            generate.confounders="sva",
                            n.confounders =10, # for the test beta matrix, I used 5 svs.
                            data=samplepheno,
                            method="rlm",
                            family="gaussian")


# show top results by p.value for each model (associations are epigenome-wide significant at p < 9.4x10-8)
top.idx <- order(ewas.ret$table$p.value)[1:10]
ewas.ret$table[top.idx,]

top.idx.adj <- order(ewas.ret.adj$table$p.value)[1:10]
ewas.ret.adj$table[top.idx.adj,]

# check if svs capture batch variable ("Sentrix_ID")
# should see significant p.values for the sv ~ Sentrix_ID
fit<- lm(ewas.ret.adj$design[, "sv1"] ~ samplepheno$Sentrix_ID)
coef(summary(fit)) 
     
     
############################################
#create outputs
############################################
ewas_res<-data.frame(
  probeID=rownames(ewas.ret$table),
  coef.none=ewas.ret$table$estimate,
  se.none=ewas.ret$table$se,
  # t.none=ewas.ret$table$t,
  p.none=ewas.ret$table$p.value,
  coef.adj=ewas.ret.adj$table$estimate,
  se.adj=ewas.ret.adj$table$se,
  # t.adj=ewas.ret.adj$table$t,
  p.adj=ewas.ret.adj$table$p.value,
  n=ewas.ret.adj$table$n)

# Adjustment for multiple testing
ewas_res$fdr.none<-p.adjust(ewas_res$p.none, method="fdr") 
ewas_res$bonferroni.none<-p.adjust(ewas_res$p.none, method="bonferroni")
ewas_res$fdr.adj<-p.adjust(ewas_res$p.adj, method="fdr") 
ewas_res$bonferroni.adj<-p.adjust(ewas_res$p.adj, method="bonferroni")

#Add information about EWAS
ewas_res$formula.none<-paste(ewas.ret$formula)
ewas_res$formula.adj<-paste(ewas.ret.adj$formula)
ewas_res$nSVs_adj<-ewas.ret$n.confounders #Number 

head(ewas_res)
 
#Add Lambda (a measure of test statistic inflation)
Lambda<-function(P){
  chisq <- qchisq(1-P,1)
  median(chisq,na.rm=T)/qchisq(0.5,1)
}

ewas_res$lambda.none<-Lambda(ewas_res$p.none)
ewas_res$lambda.adj<-Lambda(ewas_res$p.adj)


#Add annotation information about probes
annotation <- meffil.get.features("epic2")
ewas_res<-merge(ewas_res,annotation,by.x="probeID", by.y="name",all.x=TRUE)
ewas_res<-ewas_res[order(ewas_res$p.adj),] # p.value of the adjusted model w/o multiple-testing correction

saved_ewas_results <- paste("./EWAS/outputs/", EWAS_name ,Sys.Date(),".csv", sep = "_")
write.csv(ewas_res, file=saved_ewas_results)


#annotation
temp<-annotation[annotation$name %in% rownames(ewas.ret$table)==T,]
m <- match(rownames(ewas.ret$table),temp$name)
temp<- temp[m,]
identical(temp$name, rownames(ewas.ret$table))

#report unadjusted model
sum.ret <- ewaff.summary(ewas.ret, temp$chr, temp$pos, betas)

ewaff.report(sum.ret, output.file=paste("./EWAS/outputs/ewaff_rep_unadj/report_", EWAS_name, "_unadj.html", sep=""),
           author="DJQ",
            study=paste(EWAS_name, "_unadj", sep=""))

##report adjusted model
sum.ret <- ewaff.summary(ewas.ret.adj, temp$chr, temp$pos, betas)

ewaff.report(sum.ret, output.file=paste("./EWAS/outputs/ewaff_rep_adj/report_", EWAS_name, "_adj.html", sep=""),
            author="DJQ",
            study=paste(EWAS_name, "_adj", sep=""))


####END OF EWAS #####

results<- ewas_res

save(ewas_res, betas, annotation, samplepheno, file=paste("EWAS/datasets_",EWAS_name ,".RData"))

############################################
# Define plot functions
############################################
#define plot functions required (adapted from meffil)

scatter.thinning <- function(x,y,resolution=100,max.per.cell=100) {
  x.cell <- floor((resolution-1)*(x - min(x,na.rm=T))/diff(range(x,na.rm=T))) + 1
  y.cell <- floor((resolution-1)*(y - min(y,na.rm=T))/diff(range(y,na.rm=T))) + 1
  z.cell <- x.cell * resolution + y.cell
  frequency.table <- table(z.cell)
  frequency <- rep(0,max(z.cell, na.rm=T))
  frequency[as.integer(names(frequency.table))] <- frequency.table
  f.cell <- frequency[z.cell]
  
  big.cells <- length(which(frequency > max.per.cell))
  sort(c(which(f.cell <= max.per.cell),
         sample(which(f.cell > max.per.cell),
                size=big.cells * max.per.cell, replace=F)),
       decreasing=F)
}



#recode X and Y chromosomes
results$chromosome<-gsub("chrX", "chr23", results$chr)
results$chromosome<-gsub("chrY", "chr24", results$chr)


#need cumulative position of each CpG in the genome (borrow meffil code)
chromosomes <- sort(unique(as.numeric(gsub("chr", "", results$chromosome))))
chromosomes<-paste("chr", chromosomes, sep="")
chromosome.lengths <- sapply(chromosomes, function(chromosome)
  max(results$position[which(results$chromosome == chromosome)]))
chromosome.lengths <- as.numeric(chromosome.lengths)
chromosome.starts <- c(1,cumsum(chromosome.lengths)+1)
names(chromosome.starts) <- c(chromosomes, "NA")
results$globalBP <- results$position + chromosome.starts[results$chromosome] - 1


#Then we need to prepare the X axis. Instead of displaying the cumulative bp position show the chromosome name instead.
axisdf<-data.frame()

j<-0
for (i in chromosomes){
  j<-j+1
  axisdf[j,1]<-i
  temp<-results[results$chromosome==i,]
  #axisdf[j,2]<-(max(temp$globalBP) + min(temp$globalBP) )/ 2
  axisdf[j,2]<-median(temp$globalBP)
}
colnames(axisdf)<-c("chromosome", "centre")


#Exp-wide EPIC-V2 threshold
sig.threshold<-9.4e-08 # 0.05/dim(results)[1]

# how many top hits
table(results$p.none < sig.threshold)
table(results$p.adj < sig.threshold)


#filter results to only plot samples with smaller p values (quicker plotting)
#results<- filter(results, -log10(results$p.none)>1)


############################################
# Manhattan plot of EWAS -unadjusted
############################################

#apply scatter thinning
selection.idx <- scatter.thinning(results$globalBP, -log10(results$p.none),
                                  resolution=100, max.per.cell=100)


#set colours
## set different patterns for CpGs in odd vs even chromosomes, and for those significant signals that are hypo vs hypermethylated
results$flag<-"odd"
results$flag[results$chromosome=="chr2" | results$chromosome=="chr4" | results$chromosome=="chr6" |
               results$chromosome=="chr8" | results$chromosome=="chr10" | results$chromosome=="chr12" |
               results$chromosome=="chr14" | results$chromosome=="chr16" | results$chromosome=="chr18" |
               results$chromosome=="chr20" | results$chromosome=="chr22" | results$chromosome=="chr24"] <-"even"
results$flag[results$coef.none>0 & results$p.none <sig.threshold]<-"hypermethylated"
results$flag[results$coef.none<0 & results$p.none <sig.threshold]<-"hypomethylated"


#Make the plot using ggplot2:

# pdf(file= paste("./EWAS/plots/",EWAS_name, "manhattan_unadj.pdf"), width=14)
plot1<-ggplot(results[selection.idx,], aes(x=globalBP, y=-log10(p.none))) +
      # Show all points
      geom_point( aes(color=as.factor(flag)), alpha=0.8, size=1.3) +
      scale_color_manual(values = c("dark grey", "red", "blue", "light grey")) +
      # custom X axis:
      scale_x_continuous( label = axisdf$chromosome, breaks= axisdf$centre ) +
      #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
      # Custom the theme:
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_text(size=5)
      ) +
      # add significance line
      geom_hline(yintercept=-log(sig.threshold,10), colour="red") +
      #add axis labels
      xlab("position") + ylab("-log10 p-value")

# plot1

# dev.off()


############################################
# QQ plot of EWAS unadjusted
############################################

# specify functions to create plot (adapted from meffil)

qq.lambda <- function(p.values, method="median", B=100) {
  stopifnot(method %in% c("median","regression","robust"))
  p.values <- na.omit(p.values)
  observed <- qchisq(p.values, df=1, lower.tail = FALSE)
  observed <- sort(observed)
  expected <- qchisq(ppoints(length(observed)), df=1, lower.tail=FALSE)
  expected <- sort(expected)
  
  lambda <- se <- NA
  if (method == "median")  {
    lambda <- median(observed)/qchisq(0.5, df=1)
    boot.medians <- sapply(1:B, function(i) median(sample(observed, replace=T)))
    se <- sd(boot.medians/qchisq(0.5,df=1))
  } else if (method %in% c("regression","robust")) {
    if (method == "regression")
      coef.table <- summary(lm(observed ~ 0 + expected))$coeff
    else
      coef.table <- summary(rlm(observed ~ 0 + expected))$coef
    lambda <- coef.table["expected",1]
    se <- coef.table["expected", "Std. Error"]
  }
  list(method=method, estimate=lambda, se=se)
}

mod.meffil.ewas.qq.plot <- function(data,
                                    sig.threshold=9.4e-8,  # EPIC experiment-wide MT correction p-threshold
                                    sig.color="red",
                                    title="QQ plot",
                                    xlab=bquote(-log[10]("expected p-values")),
                                    ylab=bquote(-log[10]("observed p-values")),
                                    lambda.method="median") {
  
  
  #sapply(names(ewas.object$analyses), function(name) {
  p.values <- sort(data$p.none, decreasing=T)
  p.values[which(p.values < .Machine$double.xmin)] <- .Machine$double.xmin
  stats <- data.frame(is.sig=p.values < sig.threshold,
                      expected=-log(sort(ppoints(p.values),decreasing=T),10),
                      observed=-log(p.values, 10))
  lambda <- qq.lambda(p.values[which(p.values > sig.threshold)],
                      method=lambda.method)
  
  label.x <- min(stats$expected) + diff(range(stats$expected))*0.1
  label.y <- min(stats$expected) + diff(range(stats$observed))*0.9
  
  lambda.label <- paste("lambda == ", format(lambda$estimate,digits=3),
                        "%+-%", format(lambda$se, digits=3),
                        "~(", lambda.method, ")", sep="")
  
  selection.idx <- scatter.thinning(stats$observed, stats$expected,
                                    resolution=100, max.per.cell=100)
  
  lim <- range(c(0, stats$expected, stats$observed))
  sig.threshold <- format(sig.threshold, digits=3)
  
  (ggplot(stats[selection.idx,], aes(x=expected, y=observed)) + 
      geom_abline(intercept = 0, slope = 1, colour="black") +              
      geom_point(aes(colour=factor(sign(is.sig)))) +
      scale_colour_manual(values=c("black", "red"),
                          name="association",
                          breaks=c("0","1"),
                          labels=c(paste("p-value >", sig.threshold),
                                   paste("p-value <", sig.threshold))) +
      annotate(geom="text", x=label.x, y=label.y, hjust=0,
               label=lambda.label,
               parse=T) +
      xlim(lim) + ylim(lim) + 
      xlab(xlab) + ylab(ylab) +
      coord_fixed() 
    #+
    #         ggtitle(paste(title, ": ", name, sep=""))
  )
  #}, simplify=F)     
}



# generate the plot
# pdf(file= "./EWAS/plots/EWAS_qq_unadj.pdf")
plot2<-mod.meffil.ewas.qq.plot(results)
# plot2
# dev.off()


############################################
# volcano plots of EWASs unadjusted
############################################

### volcano plot analysis ###
df<-results
# df$threshold <- as.factor(abs(df$coef.none) > 0.05 & df$p.none < (0.05/length(df[,1])))
df$threshold <- as.factor(abs(df$coef.none) > 0.05 & df$p.none < sig.threshold )

selection.idx <- scatter.thinning(df$coef.none, -log10(df$p.none),
                                  resolution=100, max.per.cell=100)

# pdf(file= "./EWAS/plots/EWAS_volcano_unadj.pdf")
V_plot = ggplot(data=df[selection.idx,], aes(x=coef.none, y=-log10(p.none), colour=threshold)) +
  geom_point(alpha=0.8, size=1.75) +
  theme(legend.position = "none") + 
  theme(panel.background=element_rect(fill="white"))+ 
  theme(panel.grid.minor=element_line(colour="white"))+
  theme(axis.line=element_line(colour="black"))+
  #theme (axis.text= element_text(size=12), axis.title=element_text(size=14,face="bold"))+
  #theme(plot.title = element_text(size = 16, face="bold"))+
  #xlim(c(-0.75, 0.75)) + ylim(c(0, 30)) +
  xlab("coefficient")+
  ylab("-log10 p-value")
V_plot<-V_plot + geom_hline(aes(yintercept=8),colour="gray29", linetype="dashed")
V_plot<-V_plot+  geom_vline(aes(xintercept=-0.05),colour="gray29", linetype="dashed")
V_plot<-V_plot+  geom_vline(aes(xintercept=0.05),colour="gray29", linetype="dashed")
V_plot<- V_plot+labs(title="")
V_plot<-V_plot + scale_colour_manual(values = c("light grey","orange"))
# V_plot
# dev.off()

rm(df)
 
############################################
# plots in one figure unadjusted
############################################

pdf(file= paste("./EWAS/plots/Allplots_", EWAS_name, "_unadjusted.pdf", sep=""), width=14)

plot1<-plot1+ ggtitle("A") +
  theme(plot.title = element_text(hjust = 0))
plot2<-plot2 + ggtitle("B") +
  theme(plot.title = element_text(hjust = 0))
V_plot<-V_plot + ggtitle("C") +
  theme(plot.title = element_text(hjust = 0))


grid.arrange(plot1, arrangeGrob(plot2, V_plot, ncol=2), ncol=1)
dev.off()


######## Unadjusted PLOTS END ####################################

############################################
# Manhattan plot of EWAS adjusted
############################################
#define plot functions required (adapted from meffil)

#apply scatter thinning
selection.idx <- scatter.thinning(results$globalBP, -log10(results$p.adj),
                                  resolution=100, max.per.cell=100)


#set colours
results$flag<-"odd"
results$flag[results$chromosome=="chr2" | results$chromosome=="chr4" | results$chromosome=="chr6" |
               results$chromosome=="chr8" | results$chromosome=="chr10" | results$chromosome=="chr12" |
               results$chromosome=="chr14" | results$chromosome=="chr16" | results$chromosome=="chr18" |
               results$chromosome=="chr20" | results$chromosome=="chr22" | results$chromosome=="chr24"] <-"even"
results$flag[results$coef.adj>0 & results$p.adj <sig.threshold]<-"hypermethylated"
results$flag[results$coef.adj<0 & results$p.adj <sig.threshold]<-"hypomethylated"


#Make the plot using ggplot2:

#pdf(file= "./EWAS/plots/ewas_manhattan.pdf", width=14)
plot1<-ggplot(results[selection.idx,], aes(x=globalBP, y=-log10(p.adj))) +
  
  # Show all points
  geom_point( aes(color=as.factor(flag)), alpha=0.8, size=1.3) +
  scale_color_manual(values = c("dark grey", "red", "blue", "light grey")) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome, breaks= axisdf$centre ) +
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=5)
  ) +
  
  # add significance line
  geom_hline(yintercept=-log(sig.threshold,10), colour="red") +
  
  #add axis labels
  xlab("position") + ylab("-log10 p-value")
#plot1

#dev.off()


############################################
# QQ plot of EWAS adjusted
############################################
# specify functions to create plot (adaptedd from meffil)

qq.lambda <- function(p.values, method="median", B=100) {
  stopifnot(method %in% c("median","regression","robust"))
  p.values <- na.omit(p.values)
  observed <- qchisq(p.values, df=1, lower.tail = FALSE)
  observed <- sort(observed)
  expected <- qchisq(ppoints(length(observed)), df=1, lower.tail=FALSE)
  expected <- sort(expected)
  
  lambda <- se <- NA
  if (method == "median")  {
    lambda <- median(observed)/qchisq(0.5, df=1)
    boot.medians <- sapply(1:B, function(i) median(sample(observed, replace=T)))
    se <- sd(boot.medians/qchisq(0.5,df=1))
  } else if (method %in% c("regression","robust")) {
    if (method == "regression")
      coef.table <- summary(lm(observed ~ 0 + expected))$coeff
    else
      coef.table <- summary(rlm(observed ~ 0 + expected))$coef
    lambda <- coef.table["expected",1]
    se <- coef.table["expected", "Std. Error"]
  }
  list(method=method, estimate=lambda, se=se)
}

mod.meffil.ewas.qq.plot <- function(data,
                                    sig.threshold=9.4e-8, # EPIC experiment-wide MT correction p-threshold
                                    sig.color="red",
                                    title="QQ plot",
                                    xlab=bquote(-log[10]("expected p-values")),
                                    ylab=bquote(-log[10]("observed p-values")),
                                    lambda.method="median") {
  
  
  #sapply(names(ewas.object$analyses), function(name) {
  p.values <- sort(data$p.adj, decreasing=T)
  p.values[which(p.values < .Machine$double.xmin)] <- .Machine$double.xmin
  stats <- data.frame(is.sig=p.values < sig.threshold,
                      expected=-log(sort(ppoints(p.values),decreasing=T),10),
                      observed=-log(p.values, 10))
  lambda <- qq.lambda(p.values[which(p.values > sig.threshold)],
                      method=lambda.method)
  
  label.x <- min(stats$expected) + diff(range(stats$expected))*0.1
  label.y <- min(stats$expected) + diff(range(stats$observed))*0.9
  
  lambda.label <- paste("lambda == ", format(lambda$estimate,digits=3),
                        "%+-%", format(lambda$se, digits=3),
                        "~(", lambda.method, ")", sep="")
  
  selection.idx <- scatter.thinning(stats$observed, stats$expected,
                                    resolution=100, max.per.cell=100)
  
  lim <- range(c(0, stats$expected, stats$observed))
  sig.threshold <- format(sig.threshold, digits=3)
  
  (ggplot(stats[selection.idx,], aes(x=expected, y=observed)) + 
      geom_abline(intercept = 0, slope = 1, colour="black") +              
      geom_point(aes(colour=factor(sign(is.sig)))) +
      scale_colour_manual(values=c("black", "red"),
                          name="association",
                          breaks=c("0","1"),
                          labels=c(paste("p-value >", sig.threshold),
                                   paste("p-value <", sig.threshold))) +
      annotate(geom="text", x=label.x, y=label.y, hjust=0,
               label=lambda.label,
               parse=T) +
      xlim(lim) + ylim(lim) + 
      xlab(xlab) + ylab(ylab) +
      coord_fixed() 
    #+
    #         ggtitle(paste(title, ": ", name, sep=""))
  )
  #}, simplify=F)     
}



# generate the plot
#pdf(file= "./EWAS/plots/EWAS_qq_adj.pdf")
plot2<-mod.meffil.ewas.qq.plot(results)
#plot2
#dev.off()


############################################
# volcano plots of EWASs adjusted
############################################

### volcano plot analysis adjusted for smoking age cells and genetic PCs ###
df<-results
df$threshold <- as.factor(abs(df$coef.adj) > 0.05 & df$p.adj < sig.threshold) # (0.05/length(df[,1]))

selection.idx <- scatter.thinning(df$coef.adj, -log10(df$p.adj),
                                  resolution=100, max.per.cell=100)

#pdf(file= "./EWAS/plots/EWAS_volcano_adj.pdf")
V_plot = ggplot(data=df[selection.idx,], aes(x=coef.adj, y=-log10(p.adj), colour=threshold)) +
  geom_point(alpha=0.8, size=1.75) +
  theme(legend.position = "none") + 
  theme(panel.background=element_rect(fill="white"))+ 
  theme(panel.grid.minor=element_line(colour="white"))+
  theme(axis.line=element_line(colour="black"))+
  #theme (axis.text= element_text(size=12), axis.title=element_text(size=14,face="bold"))+
  #theme(plot.title = element_text(size = 16, face="bold"))+
  #xlim(c(-0.75, 0.75)) + ylim(c(0, 30)) +
  xlab("coefficient")+
  ylab("-log10 p-value")
V_plot<-V_plot + geom_hline(aes(yintercept=8),colour="gray29", linetype="dashed")
V_plot<-V_plot+  geom_vline(aes(xintercept=-0.05),colour="gray29", linetype="dashed")
V_plot<-V_plot+  geom_vline(aes(xintercept=0.05),colour="gray29", linetype="dashed")
V_plot<- V_plot+labs(title="")
V_plot<-V_plot + scale_colour_manual(values = c("light grey","orange"))
#V_plot
#dev.off()

rm(df) 

############################################
# plots in one figure adjusted
############################################

pdf(file= paste("./EWAS/plots/Allplots_", EWAS_name, "_adjusted", ".pdf", sep=""), width=14)

plot1<-plot1+ ggtitle("A") +
  theme(plot.title = element_text(hjust = 0))
plot2<-plot2 + ggtitle("B") +
  theme(plot.title = element_text(hjust = 0))
V_plot<-V_plot + ggtitle("C") +
  theme(plot.title = element_text(hjust = 0))


grid.arrange(plot1, arrangeGrob(plot2, V_plot, ncol=2), ncol=1)
dev.off()


######## adjusted PLOTS END ####################################



######## END ####################################
