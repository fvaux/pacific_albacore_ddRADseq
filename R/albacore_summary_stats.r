#adapted code from ELIZABETH LEE
#install.packages("remotes")
#remotes::install_github("zakrobinson/RLDNe")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SNPRelate")

###CLEAR BUFFERS
rm(list=ls())
###SET DIRECTORY
setwd("C:/Users/vauxf/Documents/R/osu-rdata")
#library(RLDNe)

library(car)
library(dartR)
library(diveRsity)
library(hierfstat)
library(pegas)
#library(rgl)
library(StAMPP)
library(qvalue)
library(vcfR)
library(vegan)

#ori-pac-1-p12-r90-a3
vcf <- vcfR::read.vcfR("ori-pac-1-p12-r90-a3.vcf")
#ori-pac-1-p12-r90-a3 12 sample areas
popInFile <- read.csv("ori-pac-1-p12-r90-a3-12P.csv", header=FALSE)
#12 sample areas
popNames <- c("BA", "CS","CN","OR","WA","BC","HW","HN","JP","PH","NC","TS")

#dn-pac-1-p12-r90-a3
vcf <- vcfR::read.vcfR("dn-pac-1-p12-r90-a3.vcf")
#ori-pac-1-p12-r90-a3 12 sample areas
popInFile <- read.csv("ori-pac-1-p12-r90-a3-12P.csv", header=FALSE)
#12 sample areas
popNames <- c("BA", "CS","CN","OR","WA","BC","HW","HN","JP","PH","NC","TS")

#ori-pac-1-p2-r95-a3
vcf <- vcfR::read.vcfR("ori-pac-1-p2-r95-a3.vcf")
#ori-pac-1-p2-r95-a3 North vs South Pacific
popInFile <- read.csv("ori-pac-1-p2-r95-a3-2P.csv", header=FALSE)
#North vs South Pacific
popNames <- c("NP", "SP")

#ori-pac-1-nc2-r95-a3
vcf <- vcfR::read.vcfR("ori-pac-1-nc2-r95-a3.vcf")
#ori-pac-1-nc2-r95-a3 North vs South Pacific
popInFile <- read.csv("ori-pac-1-nc2-r95-a3-2S.csv", header=FALSE)
#North vs South Pacific
popNames <- c("M", "F")


#MAKE GENIND OBJECT FROM VCF
x <- vcfR2genind(vcf)
x

#POP NAMES FOR GENIND
pop.names <- popInFile$V1
length(pop.names)
pop(x)= pop.names

#MAKE HIERFSTAT FROM GENIND
x2 <- genind2hierfstat(x) 

#SUMMARY STATISTICS OVERALL
basicstat <- basic.stats(x2, diploid = TRUE, digits = 2) 
names(basicstat)   

allelic<-allelic.richness(x2,diploid=TRUE)
colnames(allelic$Ar) <- popNames
obs.het<-data.frame(basicstat$Ho)
exp.het<-data.frame(basicstat$Hs)
pop.freq<-data.frame(basicstat$pop.freq)
fis<-data.frame(basicstat$Fis)
perloc<-data.frame(basicstat$perloc)
overall<-data.frame(basicstat$overall)
overall

#SUMMARY STATS PER POP
obs.het.Mean <- colMeans(obs.het, na.rm = FALSE, dims = 1)
obs.het.Mean
exp.het.Mean <- colMeans(exp.het, na.rm = FALSE, dims = 1)
exp.het.Mean
fis.Mean <- colMeans(fis, na.rm = TRUE, dims = 1)
fis.Mean
allelic.Mean <- colMeans(allelic$Ar, na.rm = FALSE, dims = 1)
allelic.Mean

#CONFIDENCE INTERVALS FOR FIS PER POP
fis.CI <- boot.ppfis(x2)
fis.CI