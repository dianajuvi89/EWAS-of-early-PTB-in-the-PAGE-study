library(data.table)
library(gdata)
library(dplyr)
library(readxl)


# set WD
setwd('~/Results/Overlap_signals_other_EWAS/')


# load data
results<-as.data.frame(fread("combined_topSites_allAnalyses_improvGeneAnnot_06.04.26.csv")) # 838  25

# remove results from chr0
results<- results %>% filter(chromosome!="chr0") # 836 (total boys=318; total girls=542; joint=24)
dim(results) # 836
table(results$analysis)
# boys girls joint 
# 294   518    24 

# select columns of interest
ewas<- results %>% select(probeID, gene.symbol, Gene.refGeneWithVer, analysis, coef.adj.bacon.x, coef.adj.bacon.y)
# create a second column including only the first gene name
ewas<- ewas %>% mutate(gene_1st=sub(";.*", "", Gene.refGeneWithVer))
ewas<- ewas %>% mutate(genesymbol_1st=sub(";.*", "", gene.symbol))


####################################
## 1) Fernando et al. 2015 BMC 
# umbilical cord blood, n=22 neonates, PTB (26-36 weeks), 450K, FDR < 0.05, n=1855 cpgs (n=196 pertaining only to PTB).  Selecting n=196 CpGs related only to PTB status and not those correlated with GA, Fem vs Male, oxytocin, etc.
####################################
fernando2015<-as.data.frame(read_excel('~/Supp_Fernando_2015/12864_2015_1915_MOESM3_ESM/Supplmentary-3a-Systematic differences (DMPs) Term vs Preterm.xlsx'))
dim(fernando2015)  # 196
table(fernando2015$`FDR adjusted  p value`< 0.05) # T 0
table(fernando2015$`t.test (pvalue)`< 0.05) # T 60
# overlap by CpG (none)
table(ewas$probeID %in% fernando2015$`Illumina probe ID`)# 0
# overlap by gene
table(ewas$Gene.refGeneWithVer %in% fernando2015$UCSC_RefGene_Name) # T 6
table(ewas$gene_1st %in% fernando2015$UCSC_RefGene_Name) # T 8
table(ewas$genesymbol_1st %in% fernando2015$UCSC_RefGene_Name) # 0
# comparison
genes_c<-left_join (ewas %>% filter(gene_1st %in% fernando2015$UCSC_RefGene_Name) %>% select(probeID, gene_1st, analysis,coef.adj.bacon.x,coef.adj.bacon.y), fernando2015 %>% select(`Illumina probe ID`, `Average beta value difference`, UCSC_RefGene_Name, `FDR adjusted  p value`), by=c("gene_1st"="UCSC_RefGene_Name"), relationship="many-to-many")

print(genes_c)

# save overlapped genes
write.csv(genes_c, file="overlapped_genes_PTB_EWAS_Fernando_2015.csv", row.names = F)



####################################
## 2) Weather Emily et al. 2022
## saliva samples, neonates, EPIC array, n=258 (155 ePTB < 34 weeks), Genome-wide P < 3.6x10-8
## Total genome-wide CpGs = 8870, but they don't provide all the list in the paper, 
## Selecting only the n=20 top-10 most p-value significant and top-10 with higher absolute effect estimate reported
####################################
emily2022<-as.data.frame(read_excel('~/Supp_Emily_2022/fcac056_supplementary_data/topsites_Emily2022_not-all-the-list.xlsx'))
dim(emily2022) # 19 
# overlap by CpG (none)
table(ewas$probeID %in% emily2022$Probe) # 0
# overlap in genes
table(ewas$Gene.refGeneWithVer %in% emily2022$Gene) # T 1
table(ewas$gene_1st %in% emily2022$Gene) # T 1
table(ewas$genesymbol_1st %in% emily2022$Gene) # 0
# comparison overlapped genes
genes_c<-left_join (ewas %>% filter(gene_1st %in% emily2022$Gene) %>% select(probeID, gene_1st, analysis,coef.adj.bacon.x,coef.adj.bacon.y), emily2022, by=c("gene_1st"="Gene"), relationship="many-to-many")
print(genes_c)
# save overlapped genes
write.csv(genes_c, file="overlapped_genes_GA_EWAS_Weather_2022.csv", row.names = F)


####################################
## 3) Hodge et al. 2024
## buccal cells ~ saliva, ONLY very preterm neonates in sample (N=538), majority Whites (51%) and Blacks (22%); exposures were GA and PMA vs DNAm
## Bonferroni GA-CpGs=2366  and PMA-CpGs=14979
####################################

# hodge2024
ga<-as.data.frame(read_excel('~/Supp_Hodge_2024/41598_2024_68071_MOESM2_ESM.xlsx', sheet=1,skip = 2))
ga$bonf<-p.adjust(ga$`p-value`, method="bonferroni")
table(ga$bonf < 0.05) # 2366
ga<- ga %>% filter(bonf < 0.05)

pma<-as.data.frame(read_excel('~/Supp_Hodge_2024/41598_2024_68071_MOESM2_ESM.xlsx', sheet=2, skip=2))
pma$bonf<-p.adjust(pma$`p-value`, method="bonferroni")
table(pma$bonf<0.05) # 14978
pma<- pma %>% filter(bonf < 0.05)

# overlap by CpG 
table(ewas$probeID %in% ga$CpG) # 35
table(ewas$probeID %in% pma$CpG) # 348

# comparison GA cpgs
ga_c<- left_join (ewas %>% filter(probeID %in% ga$CpG) %>% select(probeID, Gene.refGeneWithVer,analysis,coef.adj.bacon.x,coef.adj.bacon.y), ga %>% select(CpG, Coefficient, `Gene Name`), by=c("probeID"="CpG"))
# create direction variable
ga_c<- ga_c %>% mutate(direction=ifelse(coef.adj.bacon.x <0 & coef.adj.bacon.y <0 & Coefficient <0, "---", ifelse(coef.adj.bacon.x >0 & coef.adj.bacon.y >0 & Coefficient >0, "+++", ifelse(coef.adj.bacon.x <0 & coef.adj.bacon.y <0 & Coefficient >0, "--+", ifelse(coef.adj.bacon.x <0 & coef.adj.bacon.y >0 & Coefficient >0, "-++", ifelse(coef.adj.bacon.x >0 & coef.adj.bacon.y >0 & Coefficient <0,"++-", ifelse(coef.adj.bacon.x <0 & coef.adj.bacon.y >0 & Coefficient <0,"-+-", ifelse(coef.adj.bacon.x >0 & coef.adj.bacon.y <0 & Coefficient <0,"+--", ifelse(coef.adj.bacon.x >0 & coef.adj.bacon.y <0 & Coefficient >0,"+-+", NA)))))))))
# order columns
ga_c<- ga_c %>% select(probeID, Gene.refGeneWithVer, analysis, coef.adj.bacon.x, coef.adj.bacon.y, Coefficient,direction, `Gene Name`)
print(ga_c)
# save data
write.csv(ga_c, file="overlapped_cpgs_GA_Hodge_2024.csv", row.names = F)

# comparison PMA cpgs
pma_c<- left_join (ewas %>% filter(probeID %in% pma$CpG) %>% select(probeID, Gene.refGeneWithVer,analysis,coef.adj.bacon.x,coef.adj.bacon.y), pma %>% select(CpG, Coefficient, `Gene Name`), by=c("probeID"="CpG"))
# create direction variable
pma_c<- pma_c %>% mutate(direction=ifelse(coef.adj.bacon.x <0 & coef.adj.bacon.y <0 & Coefficient <0, "---", ifelse(coef.adj.bacon.x >0 & coef.adj.bacon.y >0 & Coefficient >0, "+++", ifelse(coef.adj.bacon.x <0 & coef.adj.bacon.y <0 & Coefficient >0, "--+", ifelse(coef.adj.bacon.x <0 & coef.adj.bacon.y >0 & Coefficient >0, "-++", ifelse(coef.adj.bacon.x >0 & coef.adj.bacon.y >0 & Coefficient <0,"++-", ifelse(coef.adj.bacon.x <0 & coef.adj.bacon.y >0 & Coefficient <0,"-+-", ifelse(coef.adj.bacon.x >0 & coef.adj.bacon.y <0 & Coefficient <0,"+--", ifelse(coef.adj.bacon.x >0 & coef.adj.bacon.y <0 & Coefficient >0,"+-+", NA)))))))))
table(pma_c$analysis,pma_c$direction, useNA="always")
# order columns
pma_c<- pma_c %>% select(probeID, Gene.refGeneWithVer, analysis, coef.adj.bacon.x, coef.adj.bacon.y, Coefficient,direction, `Gene Name`)
print(pma_c)
# save data
write.csv(pma_c, file="overlapped_cpgs_PMA_Hodge_2024.csv", row.names = F)


####################################
## 4) Winchester et al. 2022
# PTB triads (child, mom, dad; n=19) vs term triads (n=21); buccal cells DNAm; identification of DMRs in each group: female child (n=136 DMRs) and Moms (n=165 DMRs) had more similar DMRs; less DMRs seen for boys; Moms and dads distinct DMRs. DNAm analysed using MeDIPseq -> up to 90% of genome-wide DNA methylation  coverage vs ~ 1% - 2% covered by arrays. 
####################################

dmr_m<-as.data.frame(read_excel('~/Supp_Winchester_2022/41598_2022_7262_MOESM6_ESM.xlsx', skip=1))

dmr_d<-as.data.frame(read_excel('~/Supp_Winchester_2022/41598_2022_7262_MOESM7_ESM.xlsx', skip=1))

dmr_g<-as.data.frame(read_excel('~/Supp_Winchester_2022/41598_2022_7262_MOESM8_ESM.xlsx', skip=1))

dmr_b<-as.data.frame(read_excel('~/Supp_Winchester_2022/41598_2022_7262_MOESM9_ESM.xlsx', skip=1))

# merge all dmrs
dmr_a<-rbind(dmr_m, dmr_d, dmr_g, dmr_b)
rm(dmr_m, dmr_d, dmr_g, dmr_b)
dmr_a<- dmr_a %>% mutate(analysis=c(rep("mom", 165), rep("dad", 73), rep("girl", 136), rep("boy", 61))) 
dmr_a$gene_r<-sub(";.*", "", dmr_a$`Gene Annotation`)

# check overlap in genes EWAS with DMRs-PTB
table(ewas$Gene.refGeneWithVer %in% dmr_a$gene_r) # T 14
table(ewas$gene_1st %in% dmr_a$gene_r) # T 18  **** 
table(ewas$genesymbol_1st %in% dmr_a$gene_r) # T 2

# comparison of estimates for overlapped genes
dmr_ewas_c<-left_join (ewas %>% filter(gene_1st %in% dmr_a$gene_r) %>% select(probeID, gene_1st, analysis,coef.adj.bacon.x,coef.adj.bacon.y), dmr_a %>% select(`DMR Name`, `CpG #`, gene_r, maxLFC, `Gene Category`, analysis), by=c("gene_1st"="gene_r"), relationship =
                         "many-to-many")

print(dmr_ewas_c)
# save comparison
write.csv(dmr_ewas_c, file="overlapped_genes_PTB-DMRs_Winchester_2022.csv", row.names = F)



####################################
## Merid et al 2020 GA meta-analysis
## This study includes non-complicated pregnancies (N=3648); n= 8899 cpg sites identified with GA at Bonferroni significance (P<1.06x10-7); cord blood DNAm; meta-analysis
####################################
# merid2020
merid<-as.data.frame(read_excel('~/Supp_Merid_2020/GestationalageEWAS_450kmeta-analysisresult.xlsx'))

merid<- merid %>% filter(PVALUE_FE < 1.06e-07 & I_SQUARE<50)

# overlap
table(ewas$probeID %in% merid$CpGID) # T 16

# comparison table
ewas_merid<-left_join (ewas %>% filter(probeID %in% merid$CpGID) %>% select(probeID, gene_1st, analysis,coef.adj.bacon.x,coef.adj.bacon.y), merid %>% select(CpGID, CHR, UCSC_KnownGene, BETA_FE, PVALUE_FE), by=c("probeID"="CpGID"))
print(ewas_merid)

# save comparison
write.csv(ewas_merid, file="overlapped_cpgs_GA_Merid_2020.csv", row.names = F)

