###########################################
# PAGE
# ePTB and Newborn DNAm
# Study: PAGE - R21 ePTB
# Analysis only in Boys (refining adjusted model and including detailed cells) 
# Author: Diana Juvinao-Quintero
# Date: 11-20-2025
###########################################

# Tutorial of usig ewaff package for EWAS
# https://htmlpreview.github.io/?https://github.com/perishky/ewaff/blob/master/docs/tutorial.html

# load libraries
library(remotes)
library(meffil)
library(ewaff)
library(meffil)
library(ggplot2)
library(gridExtra)
library(data.table)
library(dplyr)
library(gdata)
library(matrixStats)
library(bacon)

# read in commands line arguments
args = commandArgs(trailingOnly=TRUE)

# save input arguments
EWAS_RESULTS_FILE = args[1]
SEX = args[2]
NORM_BETA_ROBJ = args[3]
PHENO_W_CELLS_FILE = args[4]

# save output arguments
EWAS_PLOTS_BACON_PLOTS_UNADJ = args[5]
EWAS_PLOTS_BACON_PLOTS_ADJ = args[6]
EWAS_BACON_RESULTS_DIR = args[7]
EWAS_RESULTS_BACON_ALLPLOTS_UNADJ_PATH = args[8]
EWAS_RESULTS_BACON_ALLPLOTS_ADJ_PATH = args[9]
EWAS_RESULTS_BACON_TOPSITES_PATH = args[10]

#set EWAS name
EWAS_name<-paste0("ePTB_", SEX, "_noNAbetas_rlm_bacon")

############################################
# load the data
############################################
ewas_res<-read.csv(paste(EWAS_RESULTS_FILE))
ewas_res<-ewas_res %>% select(-X)
############################################
# Apply Bacon to adjust inflation
############################################
# unadjusted model
bc.none<-bacon(NULL, effectsizes = ewas_res$coef.none, standarderrors= ewas_res$se.none )
estimates(bc.none)
inflation(bc.none) # 1.208401 
bias(bc.none)
# plots
pdf(file= paste(EWAS_PLOTS_BACON_PLOTS_UNADJ))
plot(bc.none, type="hist")
plot(bc.none, type="qq")
dev.off()

# adjusted model
bc.adj<-bacon(NULL, effectsizes = ewas_res$coef.adj, standarderrors= ewas_res$se.adj )
estimates(bc.adj)
inflation(bc.adj)
bias(bc.adj)
# plots
pdf(file= paste(EWAS_PLOTS_BACON_PLOTS_ADJ))
plot(bc.adj, type="hist")
plot(bc.adj, type="qq")
dev.off()

# add bacon-adjusted coefficients to results
#tstat(bc)
ewas_res$coef.adj.bacon <- es(bc.adj) #extract inflation- and bias-corrected effect-sizes
ewas_res$se.adj.bacon<- se(bc.adj)  #extract inflation- and bias-corrected standard errors
ewas_res$p.adj.bacon <- pval(bc.adj) #extract inflation- and bias-corrected P-values

ewas_res$coef.none.bacon <- es(bc.none) #extract inflation- and bias-corrected effect-sizes
ewas_res$se.none.bacon<- se(bc.none)  #extract inflation- and bias-corrected standard errors
ewas_res$p.none.bacon <- pval(bc.none) #extract inflation- and bias-corrected P-values


Lambda<-function(P){
  chisq <- qchisq(1-P,1)
  median(chisq,na.rm=T)/qchisq(0.5,1)
}

ewas_res$lambda.adj.bacon<-Lambda(ewas_res$p.adj.bacon)
ewas_res$lambda.none.bacon<-Lambda(ewas_res$p.none.bacon)

# order results by bacon.adj.p
ewas_res<-ewas_res[order(ewas_res$p.adj.bacon),] # p.value of the adjusted model w/o multiple-testing correction

# save updated results
saved_ewas_results <- paste(EWAS_BACON_RESULTS_DIR, EWAS_name ,Sys.Date(),".csv", sep = "_")
write.csv(ewas_res, file=saved_ewas_results)

# rename object witht results
results<- ewas_res

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
table(results$p.none.bacon < sig.threshold)
table(results$p.adj.bacon < sig.threshold)


#filter results to only plot samples with smaller p values (quicker plotting)
#results<- filter(results, -log10(results$p.none.bacon)>1)


############################################
# Manhattan plot of EWAS -unadjusted
############################################

#apply scatter thinning
selection.idx <- scatter.thinning(results$globalBP, -log10(results$p.none.bacon),
                                  resolution=100, max.per.cell=100)


#set colours
## set different patterns for CpGs in odd vs even chromosomes, and for those significant signals that are hypo vs hypermethylated
results$flag<-"odd"
results$flag[results$chromosome=="chr2" | results$chromosome=="chr4" | results$chromosome=="chr6" |
               results$chromosome=="chr8" | results$chromosome=="chr10" | results$chromosome=="chr12" |
               results$chromosome=="chr14" | results$chromosome=="chr16" | results$chromosome=="chr18" |
               results$chromosome=="chr20" | results$chromosome=="chr22" | results$chromosome=="chr24"] <-"even"
results$flag[results$coef.none.bacon>0 & results$p.none.bacon <sig.threshold]<-"hypermethylated"
results$flag[results$coef.none.bacon<0 & results$p.none.bacon <sig.threshold]<-"hypomethylated"


#Make the plot using ggplot2:

# pdf(file= paste("./EWAS/plots/",EWAS_name, "manhattan_unadj.pdf"), width=14)
plot1<-ggplot(results[selection.idx,], aes(x=globalBP, y=-log10(p.none.bacon))) +
  # Show all points
  geom_point( aes(color=as.factor(flag)), alpha=0.8, size=1.3) +
  scale_color_manual(values = c("dark grey", "red", "blue", "light grey")) +
  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome, breaks= axisdf$centre ) +
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="bacon.none",
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
  p.values <- sort(data$p.none.bacon, decreasing=T)
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
# df$threshold <- as.factor(abs(df$coef.none.bacon) > 0.05 & df$p.none.bacon < (0.05/length(df[,1])))
df$threshold <- as.factor(abs(df$coef.none.bacon) > 0.05 & df$p.none.bacon < sig.threshold )

selection.idx <- scatter.thinning(df$coef.none.bacon, -log10(df$p.none.bacon),
                                  resolution=100, max.per.cell=100)

# pdf(file= "./EWAS/plots/EWAS_volcano_unadj.pdf")
V_plot = ggplot(data=df[selection.idx,], aes(x=coef.none.bacon, y=-log10(p.none.bacon), colour=threshold)) +
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

pdf(file= paste(EWAS_RESULTS_BACON_ALLPLOTS_UNADJ_PATH), width=14)

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
selection.idx <- scatter.thinning(results$globalBP, -log10(results$p.adj.bacon),
                                  resolution=100, max.per.cell=100)


#set colours
results$flag<-"odd"
results$flag[results$chromosome=="chr2" | results$chromosome=="chr4" | results$chromosome=="chr6" |
               results$chromosome=="chr8" | results$chromosome=="chr10" | results$chromosome=="chr12" |
               results$chromosome=="chr14" | results$chromosome=="chr16" | results$chromosome=="chr18" |
               results$chromosome=="chr20" | results$chromosome=="chr22" | results$chromosome=="chr24"] <-"even"
results$flag[results$coef.adj.bacon>0 & results$p.adj.bacon <sig.threshold]<-"hypermethylated"
results$flag[results$coef.adj.bacon<0 & results$p.adj.bacon <sig.threshold]<-"hypomethylated"


#Make the plot using ggplot2:

#pdf(file= "./EWAS/plots/ewas_manhattan.pdf", width=14)
plot1<-ggplot(results[selection.idx,], aes(x=globalBP, y=-log10(p.adj.bacon))) +
  
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
  p.values <- sort(data$p.adj.bacon, decreasing=T)
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
df$threshold <- as.factor(abs(df$coef.adj.bacon) > 0.05 & df$p.adj.bacon < sig.threshold) # (0.05/length(df[,1]))

selection.idx <- scatter.thinning(df$coef.adj.bacon, -log10(df$p.adj.bacon),
                                  resolution=100, max.per.cell=100)

#pdf(file= "./EWAS/plots/EWAS_volcano_adj.pdf")
V_plot = ggplot(data=df[selection.idx,], aes(x=coef.adj.bacon, y=-log10(p.adj.bacon), colour=threshold)) +
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

pdf(file= paste(EWAS_RESULTS_BACON_ALLPLOTS_ADJ_PATH), width=14)

plot1<-plot1+ ggtitle("A") +
  theme(plot.title = element_text(hjust = 0))
plot2<-plot2 + ggtitle("B") +
  theme(plot.title = element_text(hjust = 0))
V_plot<-V_plot + ggtitle("C") +
  theme(plot.title = element_text(hjust = 0))


grid.arrange(plot1, arrangeGrob(plot2, V_plot, ncol=2), ncol=1)
dev.off()


######## adjusted PLOTS END ####################################


############################################
# Table of top sites 
############################################

# free wkspc
keep(results, sure=T)

results<- results %>% mutate(bonferroni.adj.bacon=p.adjust(p.adj.bacon, method="bonferroni"),
                             fdr.adj.bacon=p.adjust(p.adj.bacon, method="fdr"))

table(results$p.adj.bacon < 9.4e-08)       # ** select these sites as top sites
table(results$bonferroni.adj.bacon < 0.05) 
table(results$fdr.adj.bacon < 0.05)        

topsites<-results %>% filter(p.adj.bacon < 9.4e-08) 

# extract methylation values for topsites among cases
load(NORM_BETA_ROBJ)
pheno<-read.csv(PHENO_W_CELLS_FILE)

cases<- pheno %>% filter(eptb_bin=="ePTB case" & child_sex==paste(SEX) & !is.na(A2age) & !is.na(employed) & !is.na(prebmi)& !is.na(gravidity)) %>% select(studyid)
dim(cases) # should be 84 correct!
cases.b<-norm.beta[row.names(norm.beta) %in% topsites$probeID, colnames(norm.beta) %in% cases$studyid]
dim(cases.b) # 8135   84
cases.b.median<-data.frame(cpg=rownames(cases.b), median=rowMedians(cases.b))
cases.b.median<- cases.b.median[match(topsites$probeID, rownames(cases.b.median)),]
all(rownames(cases.b.median)==topsites$probeID) # T
identical(rownames(cases.b.median),topsites$probeID) # T

# extract methylation values for topsites among controls
control<- pheno %>% filter(eptb_bin=="control" & child_sex==paste(SEX) & !is.na(A2age) & !is.na(employed) & !is.na(prebmi)& !is.na(gravidity)) %>% select(studyid)
dim(control) # should and 55 correct!
control.b<-norm.beta[row.names(norm.beta) %in% topsites$probeID, colnames(norm.beta) %in% control$studyid]
dim(control.b) # 8135 55  
control.b.median<-data.frame(cpg=rownames(control.b), median=rowMedians(control.b))
control.b.median<- control.b.median[match(topsites$probeID, rownames(control.b.median)),]
all(rownames(control.b.median)==topsites$probeID) # T
identical(rownames(control.b.median),topsites$probeID) # T

# prepare Table data
topsites<- topsites %>% mutate(rowMed.cases=cases.b.median$median, rowMed.cont=control.b.median$median)
topsites<- topsites %>% select(probeID, chromosome:position, rowMed.cases, rowMed.cont, gene.symbol:snp.exclude, coef.adj:n,lambda.adj,coef.adj.bacon:p.adj.bacon, bonferroni.adj.bacon,lambda.adj.bacon,formula.adj)
head(topsites); dim(topsites) # 8135   22
topsites <- topsites %>% mutate(gene.symbol=gsub(";.*$", "", gene.symbol))
topsites<- topsites %>% arrange(p.adj.bacon)

# save results
write.csv(topsites, file=paste(EWAS_RESULTS_BACON_TOPSITES_PATH))



######## END ####################################