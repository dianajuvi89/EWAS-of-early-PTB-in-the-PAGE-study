library(data.table)
library(gdata)
library(cowplot)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(readr)
library(reshape)
library(viridis)
library(missMethyl)
library(ggrepel) # needed for plotting
library(BiocManager) 
library(gprofiler2) # for pathway analysis plotting
library(sesameData) # gProfiler pathway analysis
library(BiocFileCache)
library(ExperimentHub)
library(sesame)  # gProfiler pathway analysis
library(readxl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# set WD
  setwd('~/Results/EWAS_ePTB')


# load data
ewasb<-as.data.frame(fread("ePTB_boys_noNAbetas_noEosinos_allCells_rlm_bacon_2025-11-20_.csv")) 
ewasg<-as.data.frame(fread("ePTB_girls_noNAbetas_noEosinos_allCells_rlm_bacon.csv")) 
results<-as.data.frame(fread("combined_topSites_allAnalyses_improvGeneAnnot_06.04.26.csv")) 


# remove results from chr0
ewasb<- ewasb %>% filter(chromosome!="chr0") # 920950 
ewasg<- ewasg %>% filter(chromosome!="chr0") # 920950
results<- results %>% filter(chromosome!="chr0") # 836 (total boys=318; total girls=542)


#######################################################################
# missMethyl Gene ontology 
#######################################################################

# Install libraries
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("missMethyl")
library(missMethyl)
# install.packages("readxl")
# library(readxl)
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# load top sites 
topb<- results %>% filter(analysis %in% c("boys", "joint")) # 318
topg<- results %>% filter(analysis %in% c("girls", "joint")) # 542

## Note: I tried using the EPIC v2 manifest for the GO/KEGG analyisis in missMethyl and it did not work. I will use instead the Sesame R package and the "KnowYourCG" function for the EPICv2-specific functional analysis. Still running the code below for comparison.

## Output from missMethyl analysis: A data frame with a row for each GO or KEGG term and the following columns:
# Term: GO term if testing GO pathways
# Ont: ontology that the GO term belongs to if testing GO pathways.
# "BP" - biological process, "CC" - cellular component, "MF" - molecular function.
# Pathway: the KEGG pathway being tested if testing KEGG terms.
# N: number of genes in the GO or KEGG term
# DE: number of genes that are differentially methylated
# P.DE: p-value for over-representation of the GO or KEGG term term
# FDR: False discovery rate

# KEGG terms
gst.kegg.b <- gometh(sig.cpg = as.character(topb$probeID),    
                   all.cpg = unique(results$probeID), 
                   collection = "KEGG", 
                   array.type = "EPIC" 
) 


gst.kegg.g <- gometh(sig.cpg = as.character(topg$probeID),    
                     all.cpg = unique(results$probeID), 
                     collection = "KEGG", 
                     array.type = "EPIC" 
) 

# top KEGG results boys 
topGSA(gst.kegg.b)  # view top 20 pathways
min(gst.kegg.b$P.DE) # P.DE = 0.05
gst.kegg.b<- gst.kegg.b[order(gst.kegg.b$P.DE),]
head(gst.kegg.b)
dim(gst.kegg.b[gst.kegg.b$P.DE<0.05,]) # 1
dim(gst.kegg.b[gst.kegg.b$FDR<0.05,]) # 0
topKEGG.b<- topGSA(gst.kegg.b)
# save top results
write.csv(topKEGG.b, file="~/Results/Functional_Analysis/topKEGG_boys.csv", row.names = F)


# top KEGG results girls 
topGSA(gst.kegg.g)  # view top 20 pathways
min(gst.kegg.g$P.DE) # P.DE = 00.126
gst.kegg.g<- gst.kegg.g[order(gst.kegg.g$P.DE),]
head(gst.kegg.g)
dim(gst.kegg.g[gst.kegg.g$P.DE<0.05,]) # 0
dim(gst.kegg.g[gst.kegg.g$FDR<0.05,]) # 0
topKEGG.g<- topGSA(gst.kegg.g)
# save top results
write.csv(topKEGG.g, file="~/Results/Functional_Analysis/topKEGG_girls.csv", row.names = F)

## GO pathways
gst.go.b <- gometh(as.character(topb$probeID), 
                 all.cpg = unique(results$probeID), 
                 collection = "GO", array.type="EPIC")

gst.go.g <- gometh(as.character(topg$probeID), 
                   all.cpg = unique(results$probeID), 
                   collection = "GO", array.type="EPIC")

# top results boys
topGSA(gst.go.b) # view top 20 terms
gst.go.b <- gst.go.b[order(gst.go.b$P.DE),]
table(gst.go.b$FDR < 0.05) # 0
gst.go.b[gst.go.b$P.DE<0.05 & gst.go.b$ONTOLOGY =="BP",]

head(gst.go.b)
gst.go.b[gst.go.b$ONTOLOGY=="CC",][1:10,]
min(gst.go.b$P.DE)
gst.go.b[1:10,]
topGOterms.b<- topGSA(gst.go.b)
# save data
write.csv(topGOterms.b, file="~/Results/Functional_Analysis/topGO_boys.csv", row.names = F)


# top results girls
topGSA(gst.go.g) # view top 20 terms
gst.go.g <- gst.go.g[order(gst.go.g$P.DE),]
table(gst.go.g$FDR < 0.05) #0
gst.go.g[gst.go.g$P.DE<0.05 & gst.go.g$ONTOLOGY =="BP",]

head(gst.go.g)
gst.go.g[gst.go.g$ONTOLOGY=="CC",][1:10,]
min(gst.go.g$P.DE)
gst.go.g[1:10,]
topGOterms.g<- topGSA(gst.go)
# save data
write.csv(topGOterms.g, file="~/Results/Functional_Analysis/topGO_girls.csv", row.names = F)


# generate Forest Plots boys
## all processes
topKeggGO.b <- rbind(data.frame(Ontology = "KEGG",
                              Name = gst.kegg.b$Description[1:10],
                              N = gst.kegg.b$N[1:10],
                              DE = gst.kegg.b$DE[1:10],
                              P.DE = gst.kegg.b$P.DE[1:10]), 
                   data.frame(Ontology = "Biological process",
                              Name = gst.go.b$TERM[gst.go.b$ONTOLOGY=="BP"][1:10],
                              N = gst.go.b$N[gst.go.b$ONTOLOGY=="BP"][1:10],
                              DE = gst.go.b$DE[gst.go.b$ONTOLOGY=="BP"][1:10],
                              P.DE = gst.go.b$P.DE[gst.go.b$ONTOLOGY=="BP"][1:10]),
                   data.frame(Ontology = "Cellular component",
                              Name = gst.go.b$TERM[gst.go.b$ONTOLOGY=="CC"][1:10],
                              N = gst.go.b$N[gst.go.b$ONTOLOGY=="CC"][1:10],
                              DE = gst.go.b$DE[gst.go.b$ONTOLOGY=="CC"][1:10],
                              P.DE = gst.go.b$P.DE[gst.go.b$ONTOLOGY=="CC"][1:10]),
                   data.frame(Ontology = "Molecular function",
                              Name = gst.go.b$TERM[gst.go.b$ONTOLOGY=="MF"][1:10],
                              N = gst.go.b$N[gst.go.b$ONTOLOGY=="MF"][1:10],
                              DE = gst.go.b$DE[gst.go.b$ONTOLOGY=="MF"][1:10],
                              P.DE = gst.go.b$P.DE[gst.go.b$ONTOLOGY=="MF"][1:10])
)

head(topKeggGO.b)
topKeggGO.b$Name <- factor(topKeggGO.b$Name, levels = rev(topKeggGO.b$Name))

tiff("~/Results/Functional_Analysis/topGOKEGG_boys.tiff", width = 16, height = 8, units="in", pointsize=1,res=300, compression="lzw")
ggplot(topKeggGO.b, aes(y = Name, x = -log10(P.DE)))+
  geom_point()+
  geom_bar(stat='identity', width=0.1)+
  theme_classic()+
  facet_wrap(~Ontology, scales = "free_y")+
  geom_vline(xintercept = -log10(0.05), col="red",linetype=3)+
  xlab("-log10(p-value)")+
  ylab("")
dev.off()

# only KEGG and GO-BP
topKeggGO.b.bp <- rbind(data.frame(Ontology = "KEGG",
                                Name = gst.kegg.b$Description[1:10],
                                N = gst.kegg.b$N[1:10],
                                DE = gst.kegg.b$DE[1:10],
                                P.DE = gst.kegg.b$P.DE[1:10]),
                     data.frame(Ontology = "Biological process",
                                Name = gst.go.b$TERM[gst.go.b$ONTOLOGY=="BP"][1:10],
                                N = gst.go.b$N[gst.go.b$ONTOLOGY=="BP"][1:10],
                                DE = gst.go.b$DE[gst.go.b$ONTOLOGY=="BP"][1:10],
                                P.DE = gst.go.b$P.DE[gst.go.b$ONTOLOGY=="BP"][1:10])
)

topKeggGO.b.bp$Ontology<- sort(as.factor(topKeggGO.b.bp$Ontology), decreasing=TRUE)

tiff("~/Results/Functional_Analysis/topGOBP_KEGG_boys.tiff", width = 12, height = 6, units="in", pointsize=1,res=300, compression="lzw")
ggplot(topKeggGO.b.bp, aes(y = Name, x = -log10(P.DE)))+
  geom_point()+
  geom_bar(stat='identity', width=0.05)+
  theme_classic()+
  facet_wrap(~Ontology, scales = "free_y")+
  geom_vline(xintercept = -log10(0.05), col="red",linetype=3)+
  xlab("-log10(p-value)")+
  ylab("")
dev.off()



# generate Forest Plots girls
## all processes
topKeggGO.g <- rbind(data.frame(Ontology = "KEGG",
                                Name = gst.kegg.g$Description[1:10],
                                N = gst.kegg.g$N[1:10],
                                DE = gst.kegg.g$DE[1:10],
                                P.DE = gst.kegg.g$P.DE[1:10]), 
                     data.frame(Ontology = "Biological process",
                                Name = gst.go.g$TERM[gst.go.g$ONTOLOGY=="BP"][1:10],
                                N = gst.go.g$N[gst.go.g$ONTOLOGY=="BP"][1:10],
                                DE = gst.go.g$DE[gst.go.g$ONTOLOGY=="BP"][1:10],
                                P.DE = gst.go.g$P.DE[gst.go.g$ONTOLOGY=="BP"][1:10]),
                     data.frame(Ontology = "Cellular component",
                                Name = gst.go.g$TERM[gst.go.g$ONTOLOGY=="CC"][1:10],
                                N = gst.go.g$N[gst.go.g$ONTOLOGY=="CC"][1:10],
                                DE = gst.go.g$DE[gst.go.g$ONTOLOGY=="CC"][1:10],
                                P.DE = gst.go.g$P.DE[gst.go.g$ONTOLOGY=="CC"][1:10]),
                     data.frame(Ontology = "Molecular function",
                                Name = gst.go.g$TERM[gst.go.g$ONTOLOGY=="MF"][1:10],
                                N = gst.go.g$N[gst.go.g$ONTOLOGY=="MF"][1:10],
                                DE = gst.go.g$DE[gst.go.g$ONTOLOGY=="MF"][1:10],
                                P.DE = gst.go.g$P.DE[gst.go.g$ONTOLOGY=="MF"][1:10])
)

head(topKeggGO.g)
topKeggGO.g$Name <- factor(topKeggGO.g$Name, levels = rev(topKeggGO.g$Name))

tiff("~/Results/Functional_Analysis/topGOKEGG_girls.tiff", width = 16, height = 8, units="in", pointsize=1,res=300, compression="lzw")
ggplot(topKeggGO.g, aes(y = Name, x = -log10(P.DE)))+
  geom_point()+
  geom_bar(stat='identity', width=0.1)+
  theme_classic()+
  facet_wrap(~Ontology, scales = "free_y")+
  geom_vline(xintercept = -log10(0.05), col="red",linetype=3)+
  xlab("-log10(p-value)")+
  ylab("")
dev.off()

# # only KEGG and GO-BP
topKeggGO.g.bp <- rbind(data.frame(Ontology = "KEGG",
                                   Name = gst.kegg.g$Description[1:10],
                                   N = gst.kegg.g$N[1:10],
                                   DE = gst.kegg.g$DE[1:10],
                                   P.DE = gst.kegg.g$P.DE[1:10]),
                        data.frame(Ontology = "Biological process",
                                   Name = gst.go.g$TERM[gst.go.g$ONTOLOGY=="BP"][1:10],
                                   N = gst.go.g$N[gst.go.g$ONTOLOGY=="BP"][1:10],
                                   DE = gst.go.g$DE[gst.go.g$ONTOLOGY=="BP"][1:10],
                                   P.DE = gst.go.g$P.DE[gst.go.g$ONTOLOGY=="BP"][1:10])
)

topKeggGO.g.bp$Ontology<- sort(as.factor(topKeggGO.g.bp$Ontology), decreasing=TRUE)

tiff("~/Results/Functional_Analysis/topGOBP_KEGG_girls.tiff", width = 12, height = 6, units="in", pointsize=1,res=300, compression="lzw")
ggplot(topKeggGO.g.bp, aes(y = Name, x = -log10(P.DE)))+
  geom_point()+
  geom_bar(stat='identity', width=0.05)+
  theme_classic()+
  facet_wrap(~Ontology, scales = "free_y")+
  geom_vline(xintercept = -log10(0.05), col="red",linetype=3)+
  xlab("-log10(p-value)")+
  ylab("")
dev.off()



#######################################################################
# Enrichment Analysis with KnowYourCG in Sesame 
#######################################################################

# # download all databases
sesameDataCache()


# set WD
setwd("~/Results/Functional_Analysis/")

# universe
universe<- ewasb
# prepare input data
query.b<- topb    # manifest %>% filter(Name %in% top_unique$CpG) %>% select(IlmnID, Name) # original probe name
query.g<- topg 

# Results boys
#####################################

# Fisher's Exact Test for enrichment using original probe name
# Note: Do not include "EPIC" as the platform, the program will recognize it by the input data
res.b<-testEnrichment(query.b$probeID, universe=universe$probeID) 
head(res.b)
# head top findings
res.b %>% filter(overlap>10) %>% select(dbname, estimate, test, FDR, overlap) %>% head
res.b %>% filter(FDR<0.05) %>% select(dbname, estimate, test, FDR, overlap) %>% print
# any enrichment?
table(res.b$FDR<0.05) # 1
# table(res$FDR<0.1) # 
table(res.b$p.value<0.05) # 14
# TIDY UP "estimate" column for Dot plot
res.b$estimate<-round(res.b$estimate, 2)

# output expected:
#' estimate= fold change using log2(OR) test. 
#' one-sided or two-sided p value, depending how it was specified in the analysis
#' dbname= name of database
#' FDR correction
#' nQ = length of query set
#' nD = length of database set
#' overlap between nQ and nD
#' Results are sorted by FDR and the magnitude of the estimate

# save results of enrichment analysis
write.csv(res.b %>% filter(FDR<0.05) %>% arrange(FDR), file="topFeature_byFDR_KCYG_boys.csv", row.names = F)

# save plots
# enrichment all-databases barplot
tiff("enrichAllDBs_boys.tiff", width =10 , height = 6, units="in", pointsize=1,res=300, compression="lzw")
KYCG_plotEnrichAll(res.b) # barplot
dev.off()
# dot plot
tiff("Dotplot_boys.tiff", width =10 , height = 6, units="in", pointsize=1,res=300, compression="lzw")
KYCG_plotDot(res.b) # set max results to plot to print clear estmate legend , n=20
dev.off()
# volcano plot
tiff("Volcanoplot_boys.tiff", width =10 , height = 6, units="in", pointsize=1,res=300, compression="lzw")
KYCG_plotVolcano(res.b) 
dev.off() 

# Gene enrichmet analysis
res.b2<- testEnrichment(query.b$probeID,KYCG_buildGeneDBs(query.b$probeID, max_distance=10000))
dim(res.b2) # 369
res.b2[1:20, c("dbname","estimate","gene_name","FDR", "nQ", "nD", "overlap")]
table(res.b2$FDR < 0.05) # 317
# save results
write.csv(res.b2 %>% filter(FDR<0.05) %>% arrange(FDR), file="geneEnrichm_KCYG_boys.csv", row.names = F)

# Pathway enrichment
gostres.b <- gost(unique(topb$Gene.refGeneWithVer), organism = "hsapiens")
dim(gostres.b$result) # 57

## This plot was done online for better output https://biit.cs.ut.ee/gprofiler/gost
## input is the list of genes annotated to CpGs detected in main model with FDR < 0.05 (resb.2 by FDR < 0.05)
tiff("pathwayAnalysis_KYCG_boys.tiff", width =13 , height = 6, units="in", pointsize=1,res=300, compression="lzw")
gostplot(gostres.b, interactive=F)
dev.off()

# tidy up dataset
table_pathways.b<-gostres.b$result
table_pathways.b$source<-as.factor(table_pathways.b$source)
table_pathways.b<- table_pathways.b %>% arrange(p_value, source)
# save results
write.csv(table_pathways.b %>% select(-parents), file="pathwayAnal_KYCG_boys.csv", row.names = F)




# Results Girls
#####################################

# Fisher's Exact Test for enrichment using original probe name
# Note: Do not include "EPIC" as the platform, the program will recognize it by the input data
res.g<-testEnrichment(query.g$probeID, universe=universe$probeID) 
head(res.g)
# head top findings
res.g %>% filter(overlap>10) %>% select(dbname, estimate, test, FDR, overlap) %>% head
res.g %>% filter(FDR<0.05) %>% select(dbname, estimate, test, FDR, overlap) %>% print
# any enrichment?
table(res.g$FDR<0.05) # 2
# table(res$FDR<0.1) # 
table(res.g$p.value<0.05) # 27
# TIDY UP "estimate" column for Dot plot
res.g$estimate<-round(res.g$estimate, 2)

# output expected:
#' estimate= fold change using log2(OR) test. 
#' one-sided or two-sided p value, depending how it was specified in the analysis
#' dbname= name of database
#' FDR correction
#' nQ = length of query set
#' nD = length of database set
#' overlap between nQ and nD
#' Results are sorted by FDR and the magnitude of the estimate

# save results of enrichment analysis
write.csv(res.g %>% filter(FDR<0.05) %>% arrange(FDR), file="topFeature_byFDR_KCYG_girls.csv", row.names = F)

# save plots
# enrichment all-databases barplot
tiff("enrichAllDBs_girls.tiff", width =10 , height = 6, units="in", pointsize=1,res=300, compression="lzw")
KYCG_plotEnrichAll(res.g) # barplot
dev.off()
# dot plot
tiff("Dotplot_girls.tiff", width =10 , height = 6, units="in", pointsize=1,res=300, compression="lzw")
KYCG_plotDot(res.g) # set max results to plot to print clear estmate legend , n=20
dev.off()
# volcano plot
tiff("Volcanoplot_girls.tiff", width =10 , height = 6, units="in", pointsize=1,res=300, compression="lzw")
KYCG_plotVolcano(res.g) 
dev.off() 

# Gene enrichmet analysis
res.g2<- testEnrichment(query.g$probeID,KYCG_buildGeneDBs(query.g$probeID, max_distance=10000))
dim(res.g2) # 596
res.g2[1:20, c("dbname","estimate","gene_name","FDR", "nQ", "nD", "overlap")]
table(res.g2$FDR < 0.05) # 53
table(res.g2$FDR < 0.05 & res.g2$overlap >=2) # 30
# save results
write.csv(res.g2 %>% filter(FDR<0.05 & overlap >=2) %>% arrange(FDR), file="geneEnrichm_KCYG_girls.csv", row.names = F)

# Pathway enrichment
gostres.g <- gost(unique(topg$Gene.refGeneWithVer), organism = "hsapiens")
dim(gostres.g$result) # 139
table(gostres.g$result$source)
# GO:BP GO:CC   HPA MIRNA    TF    WP 
#     2     4     4    28   100     1 

## This plot was done online for better output https://biit.cs.ut.ee/gprofiler/gost
## input is the list of genes annotated to CpGs detected in main model with FDR < 0.05 (resb.2 by FDR < 0.05)
tiff("pathwayAnalysis_KYCG_girls.tiff", width =13 , height = 6, units="in", pointsize=1,res=300, compression="lzw")
gostplot(gostres.g, interactive=F)
dev.off()

# tidy up dataset
table_pathways.g<-gostres.g$result
table_pathways.g$source<-as.factor(table_pathways.g$source)
table_pathways.g<- table_pathways.g %>% arrange(p_value, source)
# save results
write.csv(table_pathways.g %>% select(-parents), file="pathwayAnal_KYCG_girls.csv", row.names = F)

