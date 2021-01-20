##############################################################################################
###                   Plot mean eigengene expression per treatment                         ###
##############################################################################################
library("WGCNA")
library("cluster")
library("affy")
library("plotrix")
library("sciplot")

### Setting folders:
dir.create("./Plots_mean_expression_eigengene_per_module")
folderIN <- "."
folderIN_data <- paste0(folderIN,"/Workspaces")
pathIN_data <- paste0(folderIN_data,"/Normalized_data.Rdata")
folderIN_functions <- paste0(folderIN,"/RFunctions")
pathIN_functions <- paste0(folderIN_functions,c("/module.functions_V1.11.R","/DAVID.functions_V2.3.R"))
folderOUT_EG_EXPRS <- paste0(folderIN,"/Plots_mean_expression_eigengene_per_module")
tag <- "Day1_FB_17_07_2015"

### Loading packages & own util functions:
Dummi <- sapply(pathIN_functions,source)

### Loading data prepared data:
load(pathIN_data)

# FEMALES:
datF  <- exprs(Exgirl)
phenF <- pData(Exgirl)
## Remove outlier samples:
INDF  <- which(!(rownames(phenF) %in% c("7032-07-118.CEL","7032-07-282.CEL")))
datF  <- datF[,INDF]
phenF <- phenF[INDF,]
# MALES:
datM  <- exprs(Exboy)
phenM <- pData(Exboy)
INDM  <- which(!(rownames(phenM) %in% c("7032-07-206.CEL","7032-07-203.CEL","7032-07-114.CEL","7032-07-194.CEL")))
datM  <- datM[,INDM]
phenM <- phenM[INDM,]
## Quick check:
all(rownames(phenF)==colnames(datF))
all(rownames(phenM)==colnames(datM))

### Load Network:
netF <- get(load("./Workspaces/netw_FEMALES_Day1_FB_17_07_2015.RData"))
netM <- get(load("./Workspaces/netw_MALES_Day1_FB_17_07_2015.RData"))


###ANOVAS WITH PLOTS###----
##Make looping plotting script----
##Also append FDR corrected P-values from Erik's ANOVA script to plots
##Consider also appending top5 GO hits - script for this at bottom but no use while DAVID is out

setwd("./Plots_mean_expression_eigengene_per_module/")

#Set consistent colors
bl<-c("#008BBA")
gr<-c("#D0C91F")
re<-c("#DC403B")

## FEMALES...
## Compute correlations with traits on basis of a one-way ANOVA:
MEF <- get.module.eg(dat=datF,net=netF)
MEF$LFL<-paste(phenF$Larval.food.level)
MEF$LFL<-factor(MEF$LFL,levels=c("Poor","Control","Rich"))
MEF$AFL<-paste(phenF$Adult.food.level)
MEF$AFL<-factor(MEF$AFL,levels=c("P","C","R"))

#Plot with adult diet on x-axis and larval diet as grouping variable
res_F <- run.aov.mod(dat=datF,net=netF,phen=phenF)
for(i in 1:3){
  pdf(paste("Female_",colnames(MEF[i]),"_",tag,".pdf"),useDingbats=FALSE)
  lineplot.CI(data=MEF,AFL,MEF[,i], group = LFL, cex = 1.5,
              xlab = "Adult diet", ylab = "Module eigengene",
              col = c(bl,gr,re), pch = c(16,16,16),bty="l", 
              main=paste("Female Module",colnames(MEF[i])))
  imp<-matrix(res_F[i,4:6])
  imp1<-cbind(colnames(res_F)[4:6],imp)
  colnames(imp1)<-c("Factor","FDR p-value")
  addtable2plot("topright",table=imp1)
  dev.off()
}


## MALES##----
##########
##########
## Compute correlations with traits on basis of a one-way ANOVA:
ME <- get.module.eg(dat=datM,net=netM)
ME$LFL<-paste(phenM$Larval.food.level)
ME$LFL<-factor(ME$LFL,levels=c("Poor","Control","Rich"))
ME$AFL<-paste(phenM$Adult.food.level)
ME$AFL<-factor(ME$AFL,levels=c("P","C","R"))

#Plot with adult diet on x-axis and larval diet as grouping variable
res_M <- run.aov.mod(dat=datM,net=netM,phen=phenM)
for(i in 1:23){
  pdf(paste("Male_",colnames(ME[i]),"_",tag,".pdf"),useDingbats=FALSE)
  lineplot.CI(data=ME,AFL,ME[,i], group = LFL, cex = 1.5,
              xlab = "Adult diet", ylab = "Module eigengene",
              col = c(bl,gr,re), pch = c(16,16,16),bty="l", 
              main=paste("Male Module",colnames(ME[i])))
  imp<-matrix(res_M[i,4:6])
  imp1<-cbind(colnames(res_M)[4:6],imp)
  colnames(imp1)<-c("Factor","FDR p-value")
  addtable2plot("topright",table=imp1)
  dev.off()
}