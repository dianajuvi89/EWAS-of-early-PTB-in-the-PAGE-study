# load libraries
library(data.table)
library(gdata)
library(cowplot)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(readr)
library(reshape)
library(viridis)
library(epitools)
library(forestplot)
library(Cairo)
library(readxl)

args = commandArgs(trailingOnly=TRUE)
EWAS = args[1]
SEX = args[2]
MANIFEST = args[3]
PVAL_THRESHOLD = as.numeric(args[4])
COEF_THRESHOLD = as.numeric(args[5])
OUTPUT_DIR = args[6]

# load in full EWAS
alltested_cgs<- read_csv(file = paste(EWAS))

# load annotation file from Illumina for EPICv2 array
featureData <- read_csv(file = paste(MANIFEST), skip = 7)

# combine the two
fData <- featureData[featureData$Name %in% alltested_cgs$probeID, ]
fData<- fData %>% distinct(Name, .keep_all=T) 

# check they have same number of CpG sites
dim(alltested_cgs) # 927761 3
dim(fData) # 927761 12

# select top CpGs based on thresholds
top_hits <- alltested_cgs %>% filter(p.adj.bacon < PVAL_THRESHOLD) %>% 
                filter(abs(coef.adj.bacon) > COEF_THRESHOLD) %>%  
                select(probeID)

# free up space
rm(featureData)
rm(alltested_cgs)

# define categories
fData <- fData %>% 
     select(Name,`450k_Enhancer`,Regulatory_Feature_Group,DNase_Hypersensitivity_NAME, Relation_to_UCSC_CpG_Island) %>%
     mutate(
        Enhancer = factor(`450k_Enhancer`, exclude = NULL),
        Promoter_Assoc = factor(ifelse(Regulatory_Feature_Group %in%  c("Promoter_Associated", "Promoter_Associated_Cell_type_specific"), 1, 0), exclude = NULL),
        DHS = factor(ifelse(DNase_Hypersensitivity_NAME!="", "TRUE",""), exclude = NULL) ,
        CpG_Island = factor(gsub(".*;", "",Relation_to_UCSC_CpG_Island), exclude = NULL)
     )

# Add these variables to the top_hits
top_hits <- left_join(top_hits, fData[, c("Name","Promoter_Assoc", "Enhancer", "DHS", "CpG_Island")], by=c("probeID"="Name"))

# function to do Chi squared test for enrichment
my.chisq.test <- function(x, y){
  
  tab_x <- table(x)
  tab_y <- table(y)
  p.1 <- prop.table(tab_x)
  p.2 <- prop.table(tab_y)
  
  ch <- chisq.test(table(x), p = p.2)
  return(c(p.1, p.2, ch$stat, P=ch$p.val))
}

## Enhancers
ch1  <- my.chisq.test(x=top_hits$Enhancer, y=fData$Enhancer)
print(ch1)
#           CpG Top Hits   |       All CpGs Tested   |            
#  Not enhancer | enhancer | Not enhancer | enhancer | X-squared | P 

## Promoters
ch2 <- my.chisq.test(x=top_hits$Promoter_Assoc, y=fData$Promoter_Assoc)
print(ch2)
#           CpG Top Hits   |       All CpGs Tested   |            
#  Not promoter | promoter | Not promoter | promoter | X-squared | P 

## DNA Hypersensitive sites (open chromatine)
ch3<- my.chisq.test(x = top_hits$DHS, y= fData$DHS)
print(ch3)
#           CpG Top Hits   |       All CpGs Tested   |            
#  Not hypersen | hypersen | Not hypersen | hypersen | X-squared | P 

## CpG islands relationship
ch4 <- my.chisq.test(x=top_hits$CpG_Island, y=fData$CpG_Island)
print(ch4)
#  CpG Top Hits & All CpGs Tested
# Island | N_Shelf | N_Shore | Open Sea | S_Shelf | S_Shore | Missing | X-squared  | P

########################################################################################
### PLOTS 
########################################################################################

top__words <- paste("Top hits for", SEX)
all_words <- "All sites tested"

d1 <- cbind.data.frame(probeset=factor(rep(c(top__words, all_words), 3), levels=c(top__words, all_words)), 
                       annot=factor(c(rep("Promoter", 2), rep("Enhancer", 2), rep("DHS", 2)),
                                    levels=c("Promoter", "Enhancer", "DHS")), 
                       prop = c(ch2[c(2,4)], ch1[c(2,4)], ch3[c(2,4)]))

p1 <- ggplot(d1, aes(x=annot, y=prop, group=probeset)) + 
  geom_bar(stat="identity", aes(fill=probeset), position = "dodge", width=0.7, 
           col = 'black', alpha =0.77) + 
  scale_fill_manual(values=c("#2b8cbe", "grey")) + 
  labs(x=NULL, y="Proportion") + 
  theme_bw() + 
  ggtitle("A")+
  theme(plot.title.position = "plot",
        plot.title = element_text(face="bold"))+
  theme(legend.title = element_blank(), 
        panel.grid.major.x = element_blank(), 
        text = element_text(family = "Gill Sans", size = 22))

temp_top_hits = table(top_hits$CpG_Island)/sum(unname(table(top_hits$CpG_Island)))
temp_all_cpgs = table(fData$CpG_Island)/sum(unname(table(fData$CpG_Island)))

d2 <- cbind.data.frame(probeset=factor(c(rep(top__words, 6),
                                         rep(all_words, 6)), 
                                       levels=c(top__words, 
                                                all_words)), 
                       annot=factor(c(names(temp_top_hits),names(temp_all_cpgs)),
                                    levels=c("Open Sea", "N_Shelf", "N_Shore", 
                                             "Island", "S_Shore", "S_Shelf"),
                                    labels=c("Open Sea", "N. Shelf", "N. Shore", 
                                             "Island", "S. Shore", "S. Shelf")), 
                       prop=c(unname(temp_top_hits), unname(temp_all_cpgs)))

rm(temp_top_hits, temp_all_cpgs)

p2 <- ggplot(d2, aes(x=annot, y=prop, group=probeset)) + 
  geom_bar(stat="identity", aes(fill=probeset), position = "dodge", width=0.7, 
           col ='black', alpha = 0.77) + 
  scale_fill_manual(values=c("#2b8cbe", "grey")) + 
  labs(x=NULL, y="Proportion") + 
  theme_bw() + 
  theme(plot.title.position = "plot",
        plot.title = element_text(face="bold"))+
  ggtitle("B") + 
  theme(legend.title = element_blank(), 
        panel.grid.major.x = element_blank(), 
        text = element_text(family = "Gill Sans", size = 22))

# save plot
tiff(paste0(OUTPUT_DIR, "Lussier_Analysis_", SEX, ".tiff"), width =14 , height = 9, units="in", pointsize=1,res=300)
grid.arrange(p1, p2, ncol = 1)
dev.off()  

