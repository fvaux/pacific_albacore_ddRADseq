###ADAPTED FROM: https://cran.r-project.org/web/packages/fsthet/vignettes/fsthet-vignette.pdf
###CREDIT: Sarah P. Flanagan (original author)

#install from github
#devtools::install_github("spflanagan/fsthet_analysis/fsthet")

#install 'fsthet' package from cran repository (NOT RECOMMENDED)
#install.packages("fsthet")

###CLEAR BUFFERS
rm(list=ls())
###SET DIRECTORY
setwd("C:/Users/vauxf/Documents/R/osu-rdata")

#load package
library(fsthet)

#import data
#Fsthet wants a standard genepop.txt (output from Stacks)

##albacore
data<-my.read.genepop("ori-pac-1-p12-r90-a3.genepop.txt")
data<-my.read.genepop("ori-pac-1-p2-r95-a3.genepop.txt")
data<-my.read.genepop("ori-pac-1-nc2-r95-a3.genepop.txt")

##~~**~~~Using the wrapper function fsthet~~~**~~##
#The above functions are all contained within the wrapper function fsthet, so you don't have to go through
#each step on its own. fsthet returns a data.frame with four columns: Locus ID, heterozygosity, Fst, and a
#True/False of whether it's an outlier
out.dat<-fsthet(data)
head(out.dat)
write.csv(out.dat,file="fsthet-output.csv",row.names=TRUE,quote=FALSE)


##~~**~~~Running components of fsthet manually~~~**~~##
#calculate actual Fst and Ht values
fsts<-calc.actual.fst(data)
head(fsts)

#Plot the actual values to see what your distribution looks like
#help("calc.actual.fst")
par(mar=c(4,4,1,1))
plot(fsts$Ht, fsts$Fst,xlab="Ht",ylab="Fst",pch=19)

##Generating quantiles
#Using boot.out
#The fst.boot function generates smoothed quantiles when you specify bootstrap=FALSE.
quant.out<-fst.boot(data, bootstrap = FALSE)

#check output
str(quant.out)
head(quant.out[[3]][[1]])

#From the results of str(quant.out), you can see that fst.boot() returns a list data.frame with three
#elements: the bootstrapped values (Fsts), the bins used in the bootstrapping (Bins), and a list of the upper
#and lower smoothed quantiles (V3).

##Plotting the results
#If you want to visualize these results, you can use plotting.cis. Plotting.cis requires the raw datapoints
#(fsts) and a list with the smoothed quantiles.
#extract the confidence interavls
quant.list<-ci.means(quant.out[[3]])
head(quant.list)

#Alternatively
#quant.list<-cis$CI0.95
#head(quant.list)

#plot the results
par(mar=c(4,4,1,1))
plotting.cis(df=fsts,ci.df=quant.list,make.file=F)

##Identifying outliers
#We can also use the find.outliers function to pull out a data.frame containing the loci that lie outside of
#the quantiles.
outliers<-find.outliers(fsts,boot.out=quant.out)
head(outliers)
write.csv(outliers,file="fsthet-outliers.csv",row.names=TRUE,quote=FALSE)

##~~**~~~Extra stuff~~~**~~##
##Look at the distribution of allele frequencies
#The analyses in fsthet use the function calc.allele.freq to calculate allele frequencies. If you're interested
#in examining the allele frequency distribution in your dataset, you can use this function on your actual data.
af.actual<-apply(data[,3:ncol(data)],2,calc.allele.freq)

#extract the minimum allele frequency for each locus
min.af<-unlist(lapply(af.actual,min))
par(mar=c(2,2,2,2))
hist(min.af)