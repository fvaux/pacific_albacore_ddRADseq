###ADAPTED FROM: https://eriqande.github.io/con-gen-2018/whoa_hands_on.nb.html
###CREDIT: Eric C. Anderson (original author)

#installing whoa
#devtools::install_github(repo = "eriqande/whoa")

###CLEAR BUFFERS
rm(list=ls())
###SET DIRECTORY
setwd("C:/Users/Felix/Documents/R/osu-rdata")

#load package
library(bindrcpp)
library(vcfR)
library(whoa)

#load data
##ALBACORE
data <- read.vcfR("ori-pac-1-p12-r90-a.vcf")
data <- read.vcfR("ori-pac-1-p12-r90-a2.vcf")
data <- read.vcfR("ori-pac-1-p12-r90-a3.vcf")

data <- read.vcfR("ori-pac-1-p2-r95-a.vcf")
data <- read.vcfR("ori-pac-1-p2-r95-a2.vcf")
data <- read.vcfR("ori-pac-1-p2-r95-a3.vcf")

data <- read.vcfR("ori-pac-1-nc2-r95-a.vcf")
data <- read.vcfR("ori-pac-1-nc2-r95-a2.vcf")
data <- read.vcfR("ori-pac-1-nc2-r95-a3.vcf")

data <- read.vcfR("dn-pac-1-p12-r90-a.vcf")

#**~~~~~~~~~WHOA ANALYSIS~~~~~~~~~**#
#first get compute expected and observed genotype frequencies
data.freqs <- exp_and_obs_geno_freqs(data)
#plot loci
#set max number of loci high enough to include all loci
geno_freqs_scatter(data.freqs, max_plot_loci = 17000)

#overall her miscall rate
data.overall <- infer_m(data, minBin = 1e15)
data.overall$m_posteriors
#5% miscall rate is good

#But note that the miscall rate may be higher at lower read depths
#First, check how many individuals and how many loci we have here:
dim(data@gt)
#left = # loci, right = # individuals (substract 1 due to header)

#FROM TUTORIAL
#e.g. So, 7382 loci and 205 individuals. That means close to 1.5 million genotypes.
#So, if we want to break that up into read depth bins, we could put 50,000 in each bin and still have a large number of bins:

#Large bins
data_binned <- infer_m(data, minBin = 100000)
posteriors_plot(data_binned$m_posteriors)

#Medium bins
data_binned <- infer_m(data, minBin = 50000)
posteriors_plot(data_binned$m_posteriors)

#Small bins
data_binned <- infer_m(data, minBin = 25000)
posteriors_plot(data_binned$m_posteriors)

#Tiny bins
data_binned <- infer_m(data, minBin = 15000)
posteriors_plot(data_binned$m_posteriors)

#**~~~~~~~~~EXTRA~~~~~~~~~**#
#Note that if you have multiple, genetically distinct populations in your VCF, 
#you can select individuals from just one population by indexing them out of the vcfR object.
#For example, if we wanted only a subset of samples from the drum VCF file we could do like this:
sams <- c("AR_001",
          "AR_003",
          "AR_004",
          "AR_005",
          "AR_008",
          "AR_010",
          "AR_012",
          "AR_014",
          "AR_015",
          "AR_016")
data_subset <- data[, c("FORMAT", sams)]
# check the dimensions of the genotypes in that object:
dim(data_subset@gt)

#Notice how you have to include the column "FORMAT" when you index the object.
#This is critical - that is the column that tells you how all the auxillary information
#that comes with the genotypes (like read depths and quality scores) is formatted.