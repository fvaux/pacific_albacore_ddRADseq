rm(list=ls())
setwd("C:/Users/vauxf/Documents/R/osu-rdata")

library(boa)
#library(coda)

# Arguments:
# - file is the name of your file ex: "output_fst.txt"
# - the q-value threshold corresponding to the target False Discovery Rate (FDR)
# - size is the size of the points and text labels for outliers
# - pos is the distance between the points and the labels 
# - highlight is a optional list of marker indices to display in red.
# - name_highlighted alows to write the indices of highlighted markers instead of using a point like the other markers
# - add_text adds the indices of the outlier markers

# Output:
# This function returns different paremeters in a list
# - outliers: the list of outliers
# - nb_outliers: the number of outliers

# Typical usage: 
# - load this file into R (file/source R code)
# - in R, go to the directory where "output_fst.txt" is (file/change current dir)
# - at the R prompt, type 

###DEFINING PLOT_BAYESCAN FUNCTION
plot_bayescan<-function(res,FDR=0.05,size=1,pos=0.35,highlight=NULL,name_highlighted=F,add_text=T)
{
  if (is.character(res))
    res=read.table(res)
  
  colfstat=5
  colq=colfstat-2
  
  highlight_rows=which(is.element(as.numeric(row.names(res)),highlight))
  non_highlight_rows=setdiff(1:nrow(res),highlight_rows)
  
  outliers=as.integer(row.names(res[res[,colq]<=FDR,]))
  
  ok_outliers=TRUE
  if (sum(res[,colq]<=FDR)==0)
    ok_outliers=FALSE;
  
  res[res[,colq]<=0.0001,colq]=0.0001
  
  # plot
  plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),xlab="log10(q value)",ylab=names(res[colfstat]),type="n")
  points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,colfstat],pch=19,cex=size)
  
  if (name_highlighted) {
    if (length(highlight_rows)>0) {
      text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
    }
  }
  else {
    points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],col="red",pch=19,cex=size)
    # add names of loci over p and vertical line
    if (ok_outliers & add_text) {
      text(log10(res[res[,colq]<=FDR,][,colq])+pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),cex=size)
    }
  }
  lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=2)
  
  return(list("outliers"=outliers,"nb_outliers"=length(outliers)))
}

###DATA FORMAT AND PREVIOUS STEPS
#1. SAVE GENEPOP FILE AS GENEPOP.TXT
#2. CONVERT GENEPOP FILE TO BAYESCAN FILE USING PGDSPIDER
#3. RUN BAYESCAN ANALYSIS (SEE COMMANDS FILE FOR RUNNING ON CGRB, OR USE GUI PROGAM OF PC)
#4. WANT .FST AND .SEL RESULTS FILES FROM BAYESCAN

###IMPORT FST DATA FROM BAYESCAN
#ALBACORE
fstdata<-read.table("ori-pac-1-p12-r90-a3-12P_fst.txt")
fstdata<-read.table("ori-pac-1-p2-r95-a3-2P_fst.txt")
fstdata<-read.table("ori-pac-1-nc2-r95-a3_fst.txt")

###PLOT BAYESCAN
# if you save the output in a variable, you can recall the different results:
#plot_bayescan(fstdata,1,FDR=0.05)
results<-plot_bayescan(fstdata,1,FDR=0.05)
results$outliers
results$nb_outliers
write.csv(results$outliers,file="outliers.csv",row.names=TRUE,quote=FALSE)


###IMPORT SEL FILE FROM BAYESCAN
# plotting posterior distribution is very easy in R with the output of BayeScan:
# first load the output file *.sel produced by BayeScan
#ALBACORE
seldata<-read.table("ori-pac-1-p12-r90-a3-12P.sel",colClasses="numeric")
seldata<-read.table("ori-pac-1-p2-r95-a3-2P.sel",colClasses="numeric")
seldata<-read.table("ori-pac-1-nc2-r95-a3.sel",colClasses="numeric")

###PLOTTING SEL FILE PARAMETERS
# choose the parameter you want to plot by setting for example:
parameter="Fst1"
# then this line will make the plot for:
plot(density(seldata[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution"))
# you can plot population specific Fst coefficient by setting
parameter="Fst1"
# if you have non-codominant data you can plot posterior for Fis coefficients in each population:
parameter="Fis1"
# if you test for selection, you can plot the posterior for alpha coefficient for selection:
parameter="alpha1"
# you also have access to the likelihood with:
parameter="logL"
# if you have the package "boa" installed, you can very easily obtain Highest Probability 
# Density Interval (HPDI) for your parameter of interest (example for the 95% interval):
boa.hpd(seldata[[parameter]],0.05)