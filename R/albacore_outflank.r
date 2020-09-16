###ADAPTED FROM: http://popgen.nescent.org/DifferentiationSNP.html
###ORIGINAL CREDIT: Stéphanie Manel (author), Zhian Kamvar (edits)
###FURTHER ADAPTED FROM NOTES BY ELIZABETH LEE

###CLEAR BUFFERS
rm(list=ls())
###SET DIRECTORY
setwd("C:/Users/vauxf/Documents/R/osu-rdata")

##################################################################
#Install and load packages
#if (!("devtools" %in% installed.packages())){install.packages(devtools)}
#library(devtools)

#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")

#if (!("vcfR" %in% installed.packages())){install.packages("vcfR")} 

#devtools::install_github("whitlock/OutFLANK")

library(OutFLANK)
library(vcfR)
library(stringr)

#######################################################################################
#****READ IN VCF****
#Code adapted from Tutorial: "Bonus: Convert VCF to OutFLANK format
#On GitHub at whilock/OutFLANK/data, you can download a vcf file of the simulations. 
#Here is a simple script to convert a vcf file into OutFLANK format, using functions 
#from the R package vcfR. Note that this code is not run with the vignette."

#Read VCF file, extract genotypes from VCF file, and change genotypes to 1,2,0 or 9 for missing data
#ALBACORE
obj.vcfR <- read.vcfR("ori-pac-1-p12-r90-a3.vcf")
obj.vcfR <- read.vcfR("ori-pac-1-p2-r95-a3.vcf")
obj.vcfR <- read.vcfR("ori-pac-1-nc2-r95-a3.vcf")

geno <- extract.gt(obj.vcfR) # Character matrix containing the genotypes
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information

G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))

G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
G[geno %in% c(NA)] <- 9

table(as.vector(G))
G

#Transpose the data to make your snp Input file
snpInFile <-t(G)

#Dimensions of the data frame should be #samples x #loci
dim(snpInFile)

#######################################################################################
#****No code changed below******

#Read in population and loci files

#loci file
#Make a two column .csv file where column one is the loci number and column two the loci ID you would like to use (can be the same as column one).
lociInFile <- read.csv("loci-ori-pac-1-p12-r90-a3.csv", header=FALSE)
head(lociInFile)

lociInFile <- read.csv("loci-ori-pac-1-p2-r95-a3.csv", header=FALSE)
head(lociInFile)

lociInFile <- read.csv("loci-ori-pac-1-nc2-r95-a3.csv", header=FALSE)
head(lociInFile)


#Dimensions of the data frame should be #loci and should match snpInFile
dim(lociInFile)
length(lociInFile$V2)
head(lociInFile$V2)

#pop file (need more than one pop)
#Make a two column .csv file where column one is each sample number and column two is the population type.
popInFile <- read.csv("pop-ori-pac-1-p12-r90-a3-12P.csv", header=FALSE)
head(popInFile)

popInFile <- read.csv("pop-ori-pac-1-p12-r90-a3-2P.csv", header=FALSE)
head(popInFile)

popInFile <- read.csv("pop-ori-pac-1-p2-r95-a3.csv", header=FALSE)
head(popInFile)

popInFile <- read.csv("pop-ori-pac-1-nc2-r95-a3.csv", header=FALSE)
head(popInFile)


#Dimensions of the data frame should be #samples and should match snpInFile
dim(popInFile)
length(popInFile$V2)
head(popInFile$V2)


#Make OutFlank file
out <- MakeDiploidFSTMat(SNPmat=snpInFile, locusNames=lociInFile$V2, popNames=popInFile$V2)

head(out)

#Plots from tutorial (https://rpubs.com/lotterhos/outflank )
par(mfrow=c(1,1))
plot(out$FST, out$FSTNoCorr, 
     xlim=c(-0.01,0.3), ylim=c(-0.01,0.3),
     pch=20)
abline(0,1)

plot(out$He, out$FSTNoCorr, pch=20, col="grey")

out1 <- OutFLANK(out, NumberOfSamples=2)

OutFLANKResultsPlotter(out1, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

hist(out1$results$pvaluesRightTail)

outlier = OutFLANK(out,NumberOfSamples = 2, 
                   RightTrimFraction = 0.05, LeftTrimFraction = 0.05,
                   qthreshold = 0.05, Hmin = 0.1)


OutFLANKResultsPlotter(outlier, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

OutFLANKResultsPlotter(outlier, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         TRUE, RightZoomFraction = 0.1, titletext = NULL)


hist(outlier$results$pvaluesRightTail)

sum(outlier$results$qvalues<0.01, na.rm=TRUE)

plot(outlier$results$He, outlier$results$FST, pch=20, col="grey")
points(outlier$results$He[outlier$results$qvalues<0.01], outlier$results$FST[outlier$results$qvalues<0.01], pch=21, col="blue")

summary(outlier$results)
head(outlier$results)

outlier$results$LocusName[outlier$results$qvalues<0.01]

write.csv(outlier$results, file = "outflank-outliers.csv")



### Note how OutFLANK identifies potential outliers at He < 0.1, even though
### these loci were excluded in the trimming algorithm
# Create logical vector for top candidates
top_candidates <- outlier$results$qvalues<0.01 & outlier$results$He>0.1

plot(outlier$results$He, outlier$results$FST, pch=20, col="grey")
points(outlier$results$He[top_candidates], outlier$results$FST[top_candidates], pch=21, col="blue")

sum(top_candidates)

# list top candidates
topcan <- outlier$results[top_candidates,]
topcan[order(topcan$LocusName),]
