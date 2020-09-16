#ORIGINAL CREDIT: ELIZABETH LEE

###CLEAR BUFFERS
rm(list=ls())
###SET WORKING DIRECTORY
setwd("C:/Users/vauxf/Documents/R/osu-rdata")

#INPUT GENEPOP AND LIST OF OUTLIER LOCI (STACKS OUTPUT GENEPOP WITH LOCI IN A ROW, TAB SEPARATED)
#STACKS2 REF MAP GENEPOP SNP POSITIONS ARE -1 FOR SOME REASON, SO EDIT OUTLIER LIST BEFOREHAND
#ALBACORE
#ori-pac-1-p12-r90-a3 (12 SAMPLE AREAS)
inputfile="ori-pac-1-p12-r90-a3.genepop.txt"
locifile="ori-pac-1-p12-r90-a3-outliers1.txt"

#ori-pac-1-p2-r95-a3 (P2 NORTH VS SOUTH PACIFIC)
inputfile="ori-pac-1-p2-r95-a3.genepop.txt"
locifile="ori-pac-1-p2-r95-a3-outliers1.txt"

#MANUALLY REMOVE 'TITLE' FROM TOP LINE OF OUTPUT FILES
#OUTPUT FILES
output="outliers_genepop.txt"
exoutput="neutral_genepop.txt" #file containing loci not in the list

#RUN
input_loci<-read.table(inputfile,skip=1,nrows=1,sep=",")
input_txt<-read.table(inputfile,skip=3,sep="\t",fill=TRUE,colClasses = "character")
genos<-input_txt[,2:ncol(input_txt)]
names(genos)<- t(input_loci)
locitomatch<-read.table(locifile)
loci<-locitomatch[,1]

newgenos<-data.frame(input_txt[1])

for (i in 1:length(loci)){
  if (loci[i] %in% names(genos)){
    newgenos<-data.frame(newgenos,genos[toString(loci[i])])
    genos[toString(loci[i])]=NULL
  }else{
    newgenos<-data.frame(newgenos,as.vector(matrix("0000",nrow=nrow(newgenos))))
  }}

writeLines("Title",output)
write.table(t(locitomatch),output,sep=",",col.names=FALSE,row.names=FALSE,quote=FALSE, append=TRUE)
write.table("pop", output,sep="\t",row.names=FALSE, quote=FALSE, col.names=FALSE, append=TRUE)
write.table(newgenos,output,sep="\t",row.names=FALSE, quote=FALSE, col.names=FALSE, append=TRUE)

exloci<-names(genos)
genos<-cbind(input_txt[1],genos)

writeLines("Title",exoutput)
write.table(t(exloci),exoutput,sep=",",row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
write.table("pop", exoutput,sep="\t",row.names=FALSE, quote=FALSE, col.names=FALSE, append=TRUE)
write.table(genos,exoutput,sep="\t",row.names=FALSE, quote=FALSE, col.names=FALSE, append=TRUE)