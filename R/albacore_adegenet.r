###ADAPTED FROM: http://popgen.nescent.org/DifferentiationSNP.html
###ORIGINAL CREDIT: Stéphanie Manel (author), Zhian Kamvar (edits)

###CLEAR BUFFERS
rm(list=ls())
###SET DIRECTORY
setwd("C:/Users/Felix/Documents/R/osu-rdata")

###LOAD LIBRARIES
library(hierfstat)
library(adegenet)
#library(ade4)
#library(pegas)
#library(constants)
library(vegan)
library(car)
library(rgl)
#library(ape)
#library(seqinr)
#library(ggplot2)

###DATA FORMAT
###PACKAGE WANTS FSTAT .DAT INPUT FILE
#1. SAVE GENEPOP FILE AS GENEPOP.TXT
#2. CONVERT GENEPOP FILE TO FSTAT .DAT FILE USING PGDSPIDER

###IMPORT DATA
##ALBACORE
####ori-pac-1-p12-r90-a3
data<-read.fstat("ori-pac-1-p12-r90-a3.dat", quiet=FALSE)

####ori-pac-1-p2-r95-a3
data<-read.fstat("ori-pac-1-p2-r95-a3.dat", quiet=FALSE)
data<-read.fstat("ori-pac-1-p2-r95-a3-ada1.dat", quiet=FALSE)
data<-read.fstat("ori-pac-1-p2-r95-a3-neu1.dat", quiet=FALSE)

####ori-pac-1-nc2-r95-a3
data<-read.fstat("ori-pac-1-nc2-r95-a3.dat", quiet=FALSE)

####dn-pac-1-p12-r90-a
data<-read.fstat("dn-pac-1-p12-r90-a.dat", quiet=FALSE)
data<-read.fstat("dn-pac-1-p12-r90-a3.dat", quiet=FALSE)

#summary stats
summary(data)
summary(seppop(data, drop = TRUE)[[3]])

data.pop = seppop(data) 
summary.by.pop = lapply(data.pop, summary) 
Hobs.ls = rep(NA, length(summary.by.pop)) 
Hexp.ls = rep(NA, length(summary.by.pop)) 

summary.by.pop$Hexp

#making a plot comparing observed heterozygosity among populations
for (i in 1:length(summary.by.pop)){ 
  Hobs.ls[i] = mean(summary.by.pop[[i]]$Hobs) 
} 
barplot(Hobs.ls, names.arg = levels(pop(data)), las = 2, main = "Observed heterozygosity", ylab = "Ho") 

data.pop <- seppop(data) 
mean.hobs <- do.call("c", lapply(data.pop, function(x) mean(summary(x)$Hobs))) 
mean.hobs[is.nan(mean.hobs)] <- NA 
barplot(mean.hobs) 

#####~~~***~~~SUMMARY STATISTICS~~~***~~~#####
#Calculate Allele Richness for FILTERED Dataset
allelic<-allelic.richness(data,diploid=TRUE)
write.table(allelic,"all.AR.txt",sep="\t")

###OBSERVED AND EXPECTED HETEROZYGOSITY: FST
# Fst following Nei (1987) on genind object
basicstats<-basic.stats(data)
basicstats
write.table(basicstats,"basicstats.txt",sep="\t")
#write.csv(basicstats,"basicstats.csv")

# Weir and Cockerham's estimate
weir<-wc(data)
weir
write.table(weir,"weir.txt",sep="\t")
#write.csv(weir,"weir.csv")

###PAIRWISE FST
weirpairwise<-genet.dist(data, method = "WC84")
weirpairwise
#write.csv(weirpairwise,"weir-pairwise.csv")

###BOOTSTRAPPED PAIRWISE FST (WHAT METHOD?) WITH CONFIDENCE INTERVALS VIA HIERFSTAT
boot.ppfst(dat=data,nboot=100,quant=c(0.025,0.975),diploid=TRUE)


#####~~~***~~~PCA FOR SNP DATA~~~***~~~#####
##ALBACORE
#PCA 12 SAMPLE AREAS
pop <- data$pop
#print(pop)
X <- tab(data, NA.method="mean")
temp <- as.integer(pop(data))
#***~~~~~~~explore only PC1 and PC2
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=2)
myCol <- transp(c("red1","red1", "red1", "goldenrod1","goldenrod1","goldenrod1","green1","green1","blue1","blue1","gray48", "gray48"),.7)[temp]
myPch <- c(19,15,17,19,15,17,19,15,19,15,19,15)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#***~~~~~~~change nf to explore additional axes (e.g. nf=3)
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=3)
myCol <- transp(c("red1","red1", "red1", "goldenrod1","goldenrod1","goldenrod1","green1","green1","blue1","blue1","gray48", "gray48"),.7)[temp]
myPch <- c(19,15,17,19,15,17,19,15,19,15,19,15)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#***~~~~~~~explore only PC1 and PC2
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=2)
myCol <- transp(c("red1","red1", "red1", "goldenrod1","goldenrod1","goldenrod1","green1","green1","blue1","blue1","black", "black"),.7)[temp]
myPch <- c(19,15,17,19,15,17,19,15,19,15,19,15)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#***~~~~~~~explore only PC1 and PC2
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=2)
myCol <- transp(c("green1","magenta", "forestgreen", "cyan","orange","lightpink","yellow","blue1","red1","gray48","black"),.7)[temp]
myPch <- c(19,15,17,19,15,17,19,15,19,15,19,15)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#myCol <- transp(c("red1","red3", "red4", "goldenrod1","goldenrod3","goldenrod4","green1","green4","blue1","deepskyblue1","gray48"),.7)[temp]
#myPch <- c(0,15,7,19,10,1,17,2,18,9,8)[temp]

#PCA 2 NORTH VS SOUTH PACIFIC
pop <- data$pop
#print(pop)
X <- tab(data, NA.method="mean")
temp <- as.integer(pop(data))
#***~~~~~~~explore only PC1 and PC2
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=2)
myCol <- transp(c("blue1","red1"),.7)[temp]
myPch <- c(19,15)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#***~~~~~~~change nf to explore additional axes (e.g. nf=3)
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=3)
myCol <- transp(c("blue1","red1"),.7)[temp]
myPch <- c(19,15)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#***~~~~~~~explore only PC1 and PC2
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=2)
myCol <- transp(c("orangered1","royalblue1"),.7)[temp]
myPch <- c(19,15)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#PCA 3 NORTH, MIXED, SOUTH ORIGIN
pop <- data$pop
#print(pop)
X <- tab(data, NA.method="mean")
temp <- as.integer(pop(data))
#***~~~~~~~explore only PC1 and PC2
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=2)
myCol <- transp(c("red1","goldenrod1","blue1"),.7)[temp]
myPch <- c(15,17,19)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#***~~~~~~~change nf to explore additional axes (e.g. nf=3)
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=3)
myCol <- transp(c("red1","goldenrod1","blue1"),.7)[temp]
myPch <- c(15,17,19)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#PCA FOR MALE VS FEMALE
pop <- data$pop
#print(pop)
X <- tab(data, NA.method="mean")
temp <- as.integer(pop(data))
#***~~~~~~~explore only PC1 and PC2
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=2)
myCol <- transp(c("dodgerblue","darkorchid1"),.7)[temp]
myPch <- c(15,19)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#***~~~~~~~change nf to explore additional axes (e.g. nf=3)
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=3)
myCol <- transp(c("dodgerblue","darkorchid1"),.7)[temp]
myPch <- c(15,19)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#PCA FOR UNKNOWN, MALE VS FEMALE
pop <- data$pop
#print(pop)
X <- tab(data, NA.method="mean")
temp <- as.integer(pop(data))
#***~~~~~~~explore only PC1 and PC2
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=2)
myCol <- transp(c("gray48","dodgerblue","darkorchid1"),.7)[temp]
myPch <- c(17,15,19)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#***~~~~~~~change nf to explore additional axes (e.g. nf=3)
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=3)
myCol <- transp(c("gray48","dodgerblue","darkorchid1"),.7)[temp]
myPch <- c(17,15,19)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#PCA FOR UNKNOWN, MALE BIAS, MALE VS FEMALE
pop <- data$pop
#print(pop)
X <- tab(data, NA.method="mean")
temp <- as.integer(pop(data))
#***~~~~~~~explore only PC1 and PC2
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=2)
myCol <- transp(c("gray48","black","dodgerblue","darkorchid1"),.7)[temp]
myPch <- c(17,15,15,19)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#***~~~~~~~change nf to explore additional axes (e.g. nf=3)
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=3)
myCol <- transp(c("gray48","black","dodgerblue","darkorchid1"),.7)[temp]
myPch <- c(17,15,15,19)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#PCA YEAR OF SAMPLING
pop <- data$pop
#print(pop)
X <- tab(data, NA.method="mean")
temp <- as.integer(pop(data))
#***~~~~~~~explore only PC1 and PC2
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=2)
myCol <- transp(c("black", "blue1", "deepskyblue1", "magenta", "forestgreen", "yellow", "orange", "red1", "darkviolet"),.7)[temp]
myPch <- c(19,19,19,19,19,19,19,19,19)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#***~~~~~~~explore only PC1 and PC2
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=2)
myCol <- transp(c("black", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#6666667"),.7)[temp]
myPch <- c(19,19,19,19,19,19,19,19,19)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#***~~~~~~~explore only PC1 and PC2
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=2)
myCol <-transp(c(rainbow(9)))
myPch <- c(19,19,19,19,19,19,19,19,19)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#***~~~~~~~explore only PC1 and PC2
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=2)
myCol <-transp(topo.colors(9))
myPch <- c(19,19,19,19,19,19,19,19,19)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)

#export PC loadings for samples
#Export loadings for multiple PCs (need to change nf= in dudi.pca function above)
write.csv(pca1$li$Axis1,file="pc1-load.csv",row.names=TRUE,quote=FALSE)
write.csv(pca1$li$Axis2,file="pc2-load.csv",row.names=TRUE,quote=FALSE)
write.csv(pca1$li$Axis3,file="pc3-load.csv",row.names=TRUE,quote=FALSE)

#Export 95% PCs for CVA
#pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE, nf=204)
#write.csv(pca1$li,file="pcs-load.csv",row.names=TRUE,quote=FALSE)

##export list of eigen values and percentage variances for each PC
#eigen values for PCs
eig.val<-pca1$eig
eig.val
#percentages of variance for each PC
eig.perc <- 100*pca1$eig/sum(pca1$eig)
eig.perc
eigen<-data.frame(eig.val,eig.perc)
eigen
#writing file with both
write.csv(eigen,file="eigen-summary.csv",row.names=TRUE,quote=FALSE)

#broken stick test to determine number of 'meaningful' PCs
#(technically no such thing as statistical significance for PCs)
#Load .csv file with listed PCs and eigen values
xx<-read.csv("eigen-summary.csv")
#xx is the eigenvalues
zz<-as.data.frame(bstick(50,tot.var=sum(xx$eig.val)))
# zz is the broken stick model
components<-seq(from=1, to=nrow(xx), by=1)
xx$comp<-components
components2<-seq(from=1, to=nrow(zz), by=1)
zz$comp<-components2
plot(xx$comp, xx$eig.val, type="h")
lines(zz$comp, zz[,1], col="red")
#If above red line = meaningful according to broken-stick test
#If below red line = non-meaningful according to broken-stick test
###Generally this corresponds to PCs that explain >5% of variance among samples

#If no PCs are meaningful, it is likely that variance among samples is low, and that there is limited genetic structuring
#There may still be weak trends or minor structure among your samples
#Important: still worth exploring non-meaningful PCs to check data, especially if there is an obvious step change in eigen values

#'allele contributions' i.e. PC loading plot for each locus
#I don't find this very useful for large SNP datasets...
#PC1
loadingplot(pca1$c1^2)

PC1_loci<-pca1$c1^2
write.csv(PC1_loci,file="PC1_loci.csv",row.names=TRUE,quote=FALSE)

#####~~~***~~~DAPC FOR SNP DATA~~~***~~~#####
###ADEGENET SERVER DAPC
adegenetServer("DAPC")

###MANUAL DAPC
##run each line separately...
#find clusters from PCs estimates from data
grp <- find.clusters(data, max.n.clust=40)
print(grp)

#DAPC
dapc1 <- dapc(data, grp$grp)
print(dapc1)
scatter(dapc1)

#optimal clusters
grp$grp
optim.a.score(dapc1)

#DAPC optimal
dapc1 <- dapc(data, grp$grp)

print(dapc1)
scatter(dapc1)

#read data
data2 <-read.genepop("alba-ori-pass7-p10-r50-a.gen")

#find number of clusters to explain data
groups <-find.clusters(data2, max.n.clust=20, choose.n.clust=FALSE, criterion='diffNgroup"')

help(find.clusters)