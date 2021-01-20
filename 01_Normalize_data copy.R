######################################################################
### Read in cel files, filter bad quality samples & normalize data ###
######################################################################
library(affy)
library(ggplot2)
library(vegan)
library(genefilter)
library(tidyverse)
library(GEOquery)

## All array data and phenotypic information is available for download at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101882
## Data can be matched to the analysis here using the "7032-07-***" identifiers in the phenotypic data

files<-c("./raw_data")
setwd(files)

##Read in CELfiles##
CELfiles<-list.celfiles()
CELfiles
pd <- read.AnnotatedDataFrame("./Master Phenotype File Both Sexes.txt",header=TRUE,row.names=1)
rawdata <- ReadAffy(filenames=CELfiles,phenoData=pd,verbose=TRUE)
setwd("..")

## Re-factor factor levels
rawdata$Time.point<-factor(rawdata$Time.point,levels=c("D.1.","10%D","90%D"))
rawdata$Adult.food.level<-factor(rawdata$Adult.food.level,levels=c("P","C","R"))
rawdata$Larval.food.level<-recode(rawdata$Larval.food.level, "P"="Poor", "C"="Control","R"="Rich")
rawdata$Sex<-factor(rawdata$Sex,levels=c("M","F"))

## Select samples from D1 of adulthood only 
rawdata <- rawdata[, rawdata$Time.point == "D.1."]

##Remove arrays that failed QC or had evidence of mating 
## Removed due to QC concerns with basic affy QC parameters
# "D.1. R/R F1": 7032-07-306
# "D.1. R/R F4": 7032-07-095
# "D.1. C/C F3" 7032-07-112

## Removed because after first module making attempt these samples had higher expression of brown module and therefore appeared mated##
# "D.1. P/C F3" : 7032-07-080
# "D.1. C/P F2" : 7032-07-122
# "D.1. P/R F1" : 7032-07-126
# "D.1. C/C F1" : 7032-07-129
# "D.1. R/P F3" : 7032-07-160

bad_array <- c("D.1. R/R F1", "D.1. R/R F4", "D.1. C/C F3", "D.1. P/C F3", "D.1. C/P F2", "D.1. P/R F1", "D.1. C/C F1", "D.1. R/P F3")
rawdata<-rawdata[, (!rawdata$Customer.Sample.Name %in% bad_array)]

## Normalize all data per sex 
Exboy<-rma(rawdata[, which(rawdata$Sex=="M")])
Exgirl<-rma(rawdata[, which(rawdata$Sex=="F")])

##Filter expression sets using these "rules" 
##remove affy quality control probes denoted by "^^AFFX")
##Do not filter out least variable genes(ie var.filter=FALSE), this would involve using completely arbitrary thresholds##
##Do not filter out probes that map to the same Entrez ID and keep only the one with most variance: keep all!
Exboy<-nsFilter(Exboy, require.entrez = FALSE, remove.dupEntrez = FALSE, var.func = IQR, 
                var.cutoff = 0.5, var.filter = FALSE, filterByQuantile=TRUE,
                feature.exclude="^AFFX")
Exboy<-Exboy$eset

Exgirl<-nsFilter(Exgirl, require.entrez = FALSE, remove.dupEntrez = FALSE, var.func = IQR, 
                 var.cutoff = 0.5, var.filter = FALSE, filterByQuantile=TRUE,
                 feature.exclude="^AFFX")
Exgirl<-Exgirl$eset

## Remove uneccesary objects and save workspace for all future work##
rm(CELfiles)
rm(bad_array)
rm(pd)
rm(rawdata)
rm(files)

save.image("./Workspaces/Normalized_data.Rdata")
