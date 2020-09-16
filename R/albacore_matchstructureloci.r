###CLEAR BUFFERS
rm(list=ls())
###SET WORKING DIRECTORY
setwd("C:/Users/vauxf/Documents/R/osu-rdata")

#!/usr/bin/env Rscript
###for importing tab-separated structure files with individuals on 2 lines###
inputfile="file.structure.tsv"
locifile="outlier_loci.txt"
output="output_structure.txt"
exoutput="otherloci_structure.txt" #file containing loci not in the list

input_loci<-read.table(inputfile,sep="\t",fill=TRUE,colClasses = "character",nrows=1)
input_txt<-read.table(inputfile,sep="\t",fill=TRUE,colClasses = "character",skip=1)
genos<-input_txt[,3:ncol(input_txt)]
names(genos)<- t(input_loci)[2:length(input_loci)]
locitomatch<-read.table(locifile)
loci<-locitomatch[,1]

newgenos<-data.frame(input_txt[1:2])

for (i in 1:length(loci)){
  if (loci[i] %in% names(genos)){
    newgenos<-data.frame(newgenos,genos[toString(loci[i])])
    genos[toString(loci[i])]=NULL
  }else{
    newgenos<-data.frame(newgenos,as.vector(matrix("0000",nrow=nrow(newgenos))))
  }}

write.table(t(locitomatch),output,sep="\t",row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(newgenos,output,sep="\t",row.names=FALSE, quote=FALSE, col.names=FALSE, append=TRUE)

exloci<-names(genos)
genos<-cbind(input_txt[1:2],genos)

write.table(t(exloci),exoutput,sep="\t",row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(genos,exoutput,sep="\t",row.names=FALSE, quote=FALSE, col.names=FALSE, append=TRUE)