#ORIGINAL CREDIT: ELIZABETH LEE

###CLEAR BUFFERS
rm(list=ls())
###SET WORKING DIRECTORY
setwd("C:/Users/vauxf/Documents/R/osu-rdata")

#INPUT STRUCTURE AND LIST OF OUTLIER LOCI (STACKS OUTPUT GENEPOP WITH INDIVIDUALS ON TWO LINES, TAB SEPARATED)
#STACKS2 REF MAP STRUCTURE SNP POSITIONS ARE -1 VS VCF, SO EDIT OUTLIER LIST BEFOREHAND
#ALBACORE
#ori-pac-1-p12-r90-a3 (12 SAMPLE AREAS)
inputfile="ori-pac-1-p12-r90-a3.structure.tsv"
locifile="ori-pac-1-p12-r90-a3-outliers1.txt"

#ori-pac-1-p2-r95-a3 (P2 NORTH VS SOUTH PACIFIC)
inputfile="ori-pac-1-p2-r95-a3.structure.tsv"
locifile="ori-pac-1-p2-r95-a3-outliers1.txt"

#OUTPUT
output="outliers_structure.txt"
exoutput="neutral_structure.txt" #file containing loci not in the list
#TWO TABS NEED TO BE ADDED BACK INTO THE HEADER ROW OF THE OUTPUT FILE (INDIVIDUAL AND POPULATION COLUMNS BEFORE GENOTYPES)

#RUN
input_loci<-read.table(inputfile,sep="\t",fill=TRUE,colClasses = "character",nrows=1)
input_txt<-read.table(inputfile,sep="\t",fill=TRUE,colClasses = "character",skip=1)
genos<-input_txt[,3:ncol(input_txt)]
names(genos)<- t(input_loci)[3:length(input_loci)]
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