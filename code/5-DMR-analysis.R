# load libraries
library(devtools)
install_github("perishky/dmrff")
library(dmrff)
library(meffil)

# read in command line arguments
args = commandArgs(trailingOnly=TRUE)
EWAS = args[1]
NORM_BETA_FILE = args[2]
SEX = args[3]

# load data
result <- read.csv(file = paste(EWAS))
annotation <- meffil.get.features("epic2")
load(paste0(NORM_BETA_FILE))

# make sure methylation data compatible with EWAS
betas <- norm.beta
rm(norm.beta)
methylation <- betas[match(result$probeID, rownames(betas)), ]

# run dmrff
dmrs <- dmrff(estimate=result$coef.adj.bacon,
              se=result$se.adj.bacon,
              p.value=result$p.adj.bacon,
              methylation=methylation,
              chr=result$chromosome,
              pos=result$position,
              maxgap=500,
              verbose=T)

# keep regions with >1 CgP site and bonferroni adjusted p<0.05
dmrs <- dmrs[which(dmrs$p.adjust < 0.05 & dmrs$n > 1),]

#find out how many CpGs in the DMRs -extract them
out<-annotation[0,]

for (i in 0:nrow(bumps.ret)-1){
  i<-i+1
  test<-annotation[as.factor(annotation$chr) %in% bumps.ret$chr[i]==T & annotation$position <=bumps.ret$end[i]
                   & annotation$position >=bumps.ret$start[i],]
  test$chr<-bumps.ret$chr[i]
  test$start<-bumps.ret$start[i]
  test$end<-bumps.ret$end[i]
  
  out<-rbind(out, test)
}

dim(out)

sum(bumps.ret$n) #number of bumps

#add info to table and save
result$in_dmr_dmrff<-"no"
result$in_dmr_dmrff[result$probe %in% out$name==T]<-"yes"

write.csv(result, "dmr_annotated.csv")

#how many in dmr and bonferroni

length(result$p.none[result$p.none<0.05/length(result$p.none) & result$in_dmr_dmrff=="yes"]) #number of Cpgs

#how many of the significant EWAS CpGs are part of a DMR:
#subset "out" to include only those where CpG are also identified in EWAS
overlap<-out[out$name %in% result$probe[result$p.none<0.05/length(result$p.none)]==T,]

#get unique DMRs that are represented by an associated EWAS CpG 
overlap<-overlap[duplicated(overlap[,c("chr", "start", "end")])==F,]

#get DMRs that are novel (ie not represented by an EWAS CpG)
test<-merge(bumps.ret, overlap, by=c("chr", "start", "end"), all.x=T)
test<-test[is.na(test$type)==T,]
dim(test)

write.csv(test, "novel_DMRs.csv")
