#CREDIT: Modified from original script by Elizabeth Lee

#Extract Depth coverage from vcf file
#Calculate and plot mean or SD for each locus

###CLEAR BUFFERS
rm(list=ls())
###SET DIRECTORY
setwd("C:/Users/Felix/Documents/R/osu-rdata")

library(vcfR)
library(matrixStats)

#Import vcf file
vcf <- vcfR::read.vcfR("FILENAME.vcf")
vcf

#Extract a matric of Depth Coverage for each sample and each locus
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
head(dp)
dim(dp)
#Save coverage depth per individual per locus to file
write.csv(dp, file = "dp.csv")

#Calculate Mean Depth coverage for each locus
RowM <- rowMeans(dp[, -ncol(dp)],na.rm = TRUE)
RowM
mean(RowM)
range(RowM)
plot(RowM)

#Calculate SD of Depth Coverage for each locus
RowSDlist <- rowSds(dp[, -ncol(dp)], na.rm=TRUE)
RowSD <- cbind(RowSDlist, dp[,0])
RowSD
mean(RowSD)
range(RowSD)
plot(RowSD)

#Make histograms and Boxplots of Depth Coverage measurements
par(mfrow=c(2,2)) 
hist(RowM, breaks=100, xlab ="Mean Coverage of Each Loci", main = "Histogram of Mean Coverage \n of Each Loci")
boxplot(RowM,main ="Mean Coverage of Loci")
hist(RowSD, breaks=100, xlab ="SD of Coverage of Each Loci", main = "Histogram of SD of Coverage \n of Each Loci")
boxplot(RowSD,main ="SD of Coverage of Loci")
par(mfrow=c(1,1)) 

#Print the outliers
#Boxplot Outliers: outside 1.5 times the interquartile range above the upper quartile and below the lower quartile
outliersM = boxplot(RowM, plot=FALSE)$out
sort(outliersM)

outliersSD = boxplot(RowSD, plot=FALSE)$out
sort(outliersSD)

#Write Means and SD to a new file
write.csv(RowSD, file = "SD_both_pops.csv")
write.csv(RowM, file = "Mean_both_pops.csv")
