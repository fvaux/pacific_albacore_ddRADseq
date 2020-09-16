###CLEAR BUFFERS
rm(list=ls())
###SET WORKING DIRECTORY
setwd("C:/Users/vauxf/Documents/R/osu-rdata")

###LOAD LIBRARIES
library(vcfR)
library(StAMPP)

###DATA FORMAT
#1. IMPORT VCF OUTPUT FILE FROM STACKS
#2. NEED TO MAKE VECTOR WITH MATCHING POPULATION NAMES, AS VCF (AND GENLIGHT OBJECTS) FILES DO NOT CONTAIN POP INFO

###IMPORT DATA
##ALBACORE
vcf<-read.vcfR("ori-pac-1-p12-r90-a3.vcf", verbose = FALSE)
vcf<-read.vcfR("ori-pac-1-p12-r90-a3-neu.vcf", verbose = FALSE)
vcf<-read.vcfR("ori-pac-1-p12-r90-a3-ada.vcf", verbose = FALSE)

vcf<-read.vcfR("ori-pac-1-p2-r95-a3.vcf", verbose = FALSE)
vcf<-read.vcfR("ori-pac-1-p2-r95-a3-neu.vcf", verbose = FALSE)
vcf<-read.vcfR("ori-pac-1-p2-r95-a3-ada.vcf", verbose = FALSE)

vcf<-read.vcfR("ori-pac-1-nc2-r95-a3.vcf", verbose = FALSE)

###MAKE POPULATION LABEL VECTORS (SAME LENGTH AS # OF INDIVIDUALS)
##ALBACORE
#308 fish, 12 sample areas
pop.names <- c("BA", "BA", "BA", "BA", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CN", "CN", "CN", "CN", "CN", "CN", "CN", "CN", "CN", "CN", "CN", "CN", "CN", "CN", "CN", "CN", "CN", "CN", "CN", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "OR", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "WA", "BC", "BC", "BC", "BC", "BC", "BC", "BC", "BC", "BC", "BC", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HW", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "HN", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "JP", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "TS", "TS", "TS", "TS", "TS", "TS", "TS", "TS", "TS", "TS", "TS", "TS", "TS", "TS", "TS", "TS", "TS", "TS", "TS", "TS")

#308 fish, North vs South Pacific
pop.names <- c("NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "NP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP", "SP")

#54 New Caledonia fish, Male vs Female
pop.names <- c("M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F")

###CONVERTING FILES
#Convert vcf into genlight object 
x <- vcfR2genlight(vcf)

#Convert genlight object to a matrix in stampp format 
x2 <- as.matrix(x) #convert genlight object to matrix
sample <- row.names(x2) #sample names

ploidy <- ploidy(x) #extract ploidy info from genlight object
x2 = x2 * (1/ploidy) #convert allele counts to frequency
x2[is.na(x2)] = NaN
format <- vector(length = length(sample))
#Format id for the genotype data
format[1:length(format)] = "freq"

x.stampp <- as.data.frame(cbind(sample, pop.names, ploidy, format, x2)) #convert to basic r data.frame suitable to stamppConvert

geno <- stamppConvert(x.stampp, 'r')
#You should now be able to run all the stampp commands with 'geno' as the input object 
#e.g. fst <- stamppFst(geno)

###*FINALLY* RUNNING STAMPP
##CALCULATING PAIRWISE FST WITH P-VALUES
data.fst<-stamppFst(geno, nboots = 10000, percent = 95, nclusters = 1)
options(max.print=40000)
tail(data.fst, 1)
head(data.fst, 2)
#names(data.fst)

#save all values
#(listed by column, not in a table)
write.csv(data.fst$Bootstraps, file = "fst-full.csv")

#options(max.print=100)
#print(data.fst)

##CALCULATE GENETIC DISTANCE BETWEEN POPULATIONS
data.D.pop<-stamppNeisD(geno, TRUE)
print(data.D.pop)

##CALCULATE GENETIC DISTANCE BETWEEN INDIVIDUALS
data.D.ind<-stamppNeisD(geno, FALSE)
#print(data.D.ind)

#CALCULATE AMOVA
amova.ind<-stamppAmova(data.D.ind, geno, 100)
print(amova.ind)

#Export the genetic distance matrix in Phylip format
#stamppPhylip(data.D.pop, file="geno-pop-d.txt")
#stamppPhylip(data.D.ind, file="geno-ind-d.txt")

#CALCULATE 'GENOMIC RELATIONSHIP VALUES' BETWEEN EACH INDIVIDUAL
#data.gmat<-stamppGmatrix(geno)
#print(data.gmat)


#############################################################################################################
#HELP FROM STAMMP CREATOR LUKE PEMBLETON
#############################################################################################################
## as per your post 
##2.1 Read vcf file into into R 
vcf<-read.vcfR("g002-p3-r60-m10.vcf", verbose = FALSE) 
##2.2 Convert vcf into genlight object 
x <- vcfR2genlight(vcf) 

### convert genlight object to a matrix in stampp format 
x2 <- as.matrix(x) #convert genlight object to matrix 
sample <- row.names(x2) #sample names 

#pop.names <- pop(x) #extract ploidy info from genlight object (however this is not available when imported via a vcf file) 
#pop.names <- #provide a vector here of the same length of the number of samples, with a corresponding population name/id
#pop.names <- c(")
#e.g. c("popA", "popA", "popB", "popB", "popB", "popC", "popC") 

ploidy <- ploidy(x) #extract ploidy info from genlight object 
x2 = x2 * (1/ploidy) #convert allele counts to frequency 
x2[is.na(x2)] = NaN 
format <- vector(length = length(sample))
#format id for the genotype data
format[1:length(format)] = "freq"  

x.stampp <- as.data.frame(cbind(sample, pop.names, ploidy, format, x2)) #convert to basic r data.frame suitable to stamppConvert 

geno <- stamppConvert(x.stampp, 'r') 

#you should now be able to run all the stampp commands with 'geno' as the input object 
# e.g. fst <- stamppFst(geno) 