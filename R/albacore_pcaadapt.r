###ADAPTED FROM: http://popgen.nescent.org/2016-01-26-SNP-selection.html
###CREDIT: Stéphanie Manel (author), Alicia Dalongeville (author), Zhian Kamvar (reviewer)

#INSTALLING PACKAGES
####PCAdapt from: https://github.com/bcm-uga/pcadapt
#install.packages("devtools")
#devtools::install_github("bcm-uga/pcadapt")
####qvalue from: https://www.bioconductor.org/packages/release/bioc/html/qvalue.html
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")

#########################********************************************************#########################
###CLEAR BUFFERS
rm(list=ls())

###SET DIRECTORY
setwd("C:/Users/vauxf/Documents/R/osu-rdata")

###LOAD LIBRARIES
library(pcadapt)
library(qvalue)

###DATA FORMAT
##PACKAGE WANTS .PED INPUT FILE
#1. CONVERT STRUCTURE.TSV TO PEDIGREE .PED FORMAT USING PGDSPIDER
#2. MANUALLY CHANGE 1, 2, 3, 4 TO A, T, C, G (0 stays 0) IN NOTEPAD++ (SAVED AS a3 OR d3.ped)

###IMPORT DATA
##ALBACORE
#308 samples (ori-pac-1-p12-r90-a4)
data<-read.pcadapt("ori-pac-1-p12-r90-a4.ped", type = "ped")

#308 samples (ori-pac-1-p2-r95-a4)
data<-read.pcadapt("ori-pac-1-p2-r95-a4.ped", type = "ped")

#54 samples (ori-pac-1-nc2-r95-a4)
data<-read.pcadapt("ori-pac-1-nc2-r95-a4.ped", type = "ped")

###PC SCREE PLOT
#25 PCs used as arbitrary number?
K <- 10
x <- pcadapt(data, K = K)
###SCREE PLOT FOR ESTIMATING K (NUMBER OF PCS OF INTEREST)
plot(x, option = "screeplot")

##REPEAT PCA WITH DECIDED K VALUE
K <- 2
x <- pcadapt(data, K = K)

###MAKING POPULATION NAMES FOR PCA
#North vs South Pacific
poplist.names <- c(rep("NP", 234),rep("SP", 74))
print(poplist.names)

#12 Pacific areas
poplist.names <- c(rep("BA", 4),rep("CS", 39),rep("CN", 19),rep("OR", 31),rep("WA", 27),rep("BC", 10),rep("HW", 25),rep("HN", 28),rep("JP", 22),rep("PH", 29),rep("NC", 54),rep("TS", 20))
print(poplist.names)

#10 North Pacific areas
poplist.names <- c(rep("BA", 4),rep("CS", 30),rep("CN", 19),rep("OR", 25),rep("WA", 22),rep("BC", 10),rep("HW", 21),rep("HN", 26),rep("JP", 22),rep("PH", 25))
print(poplist.names)

#6 West Coast areas
poplist.names <- c(rep("BA", 4),rep("CS", 30),rep("CN", 19),rep("OR", 25),rep("WA", 22),rep("BC", 10))
print(poplist.names)

#2 West Coast regions
poplist.names <- c(rep("Southern", 57),rep("Northern", 53))
print(poplist.names)

#New Caledonia males vs females
poplist.names <- c(rep("Male", 36),rep("Female", 18))
print(poplist.names)

###PCA
##without population colouring
#plot(x, option = "scores")
##with population colouring
#plot(x, option = "scores", pop = data[, 1])
plot(x, option = "scores", pop = poplist.names)

###EXAMINE FURTHER PCS (E.G. PCs 3 and 4)
plot(x, option = "scores", i = 3, j = 2, pop = poplist.names)

###EXAMINE FURTHER PCS (E.G. PCs 5 and 6)
plot(x, option = "scores", i = 1, j = 4, pop = poplist.names)


###NUMERICAL QUANTITIES OBTAINED AFTER PERFORMING A PCA
summary(x)

###MANHATTAN PLOT
plot(x, option = "manhattan")

###Q PLOT WITH PCADAPT DEFAULT THRESHOLD
plot(x, option = "qqplot", threshold = 0.1)

###Q PLOT WITH PCADAPT MORE STRINGENT THRESHOLD
plot(x, option = "qqplot", threshold = 0.05)

###Q PLOT WITH PCADAPT EVEN MORE STRINGENT THRESHOLD
plot(x, option = "qqplot", threshold = 0.01)

###DISTRIBUTION OF MAHALANOBIS DISTANCES
plot(x, option = "stat.distribution")

###ESTIMATING Q VALUES
qval <- qvalue(x$pvalues)$qvalues

###PCADAPT DEFAULT SETTINGS FOR IDENTIFYING OUTLIERS
alpha <- 0.1
outliers_pcadapt <- which(qval < alpha)
print(outliers_pcadapt)
length(outliers_pcadapt)

##MORE STRINGENT THRESHOLD TO DETECT OUTLIERS
alpha <- 0.05 
outliers <- which(qval < alpha)
print(outliers)
length(outliers)

##EVEN MORE STRINGENT THRESHOLD TO DETECT OUTLIERS
alpha <- 0.01 
outliers <- which(qval < alpha)
print(outliers)
length(outliers)