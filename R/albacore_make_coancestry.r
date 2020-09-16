###CLEAR BUFFERS
rm(list=ls())
###SET DIRECTORY
setwd("C:/Users/vauxf/Documents/R/osu-rdata")

#load packages
library(vcfR)
library(pegas)

#Write a genotype input file from a VCF file
#Read VCF file and change to a loci object type for pegas
my_vcf <- read.vcfR("FILENAME.vcf")
my_vcf
x <- vcfR2loci(my_vcf) #, return.alleles = TRUE)
x

#Write a Coancestory genotype input file called "final_genotypefile.txt"
write.loci(x, loci.sep = "\t", allele.sep = "\t", quote = FALSE, col.names = FALSE, na = "NA\tNA", file = "genotypefile.txt")
df <- read.table("genotypefile.txt", sep = "\t", header = FALSE, row.names = 1)
df[] <- lapply(df, function(x) as.integer(x)+1)
write.table(df, file = "final_genotypefile.txt", sep = "\t", na = "0", row.names=TRUE, col.names = FALSE)

#File output the working directory
########
#Make Allele Frequency files by extracting AD from VCF files on the command line with:
vcftools --vcf populations.snps.vcf --freq --out output
#Or can just use the "Calculate Allele Frequecny" setting in Coancestry

###########################
##Typical kinships to use:
.0 .0 .0 .0 .0 .0 .0 .5 .5 100 
.0 .0 .0 .0 .0 .0 .0 1. 0 100 
.0 .0 .0 .0 .0 .0 .25 .5 .25 100 
.0 .0 .0 .0 .0 .0 .25 0 .75 100 
.0 .0 .0 .0 .0 .0 .0 .25 .75 100