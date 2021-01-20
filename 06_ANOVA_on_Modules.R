##################################################################
### 14_08_2015_Rerun_analysis_from_beginning: ANOVA on MODULES ###
##################################################################
library("WGCNA")
library("cluster")
library("affy")
library("knitr")

### Setting folders:
folderIN <- "."
folderIN_data <- paste0(folderIN,"/Workspaces")
pathIN_data <- paste0(folderIN_data,"/Normalized_data.Rdata")
folderIN_functions <- paste0(folderIN,"/RFunctions")
pathIN_functions <- paste0(folderIN_functions,c("/module.functions_V1.11.R","/DAVID.functions_V2.3.R"))
folderOUT_QC <- paste0(folderIN,"/QC")
folderOUT_TISSUE <- paste0(folderIN,"/Tissue_output")
tag <- "Day1_FB_17_07_2015"

### Loading own util functions:
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

#########################
### ANOVA on MODULES ###
#########################
# FEMALES:
res_F <- run.aov.mod(dat=datF,net=netF,phen=phenF)
kable(format(res_F,digits=2))
#|            |Larval  |Adult   |LarvalxAdult |fdr_Larval |fdr_Adult |fdr_LarvalxAdult |sign_Larval |sign_Adult |sign_LarvalxAdult |
#|:-----------|:-------|:-------|:------------|:----------|:---------|:----------------|:-----------|:----------|:-----------------|
#|MEblue      |1.4e-09 |7.8e-03 |3.3e-01      |1.2e-08    |2.4e-02   |4.9e-01          |1.0e+00     |0.0e+00    |0.0e+00           | L
#|MEturquoise |2.0e-08 |1.8e-01 |5.6e-01      |9.1e-08    |3.2e-01   |6.8e-01          |1.0e+00     |0.0e+00    |0.0e+00           | L
#|MEgrey      |6.9e-01 |1.6e-02 |6.0e-01      |6.9e-01    |3.5e-02   |6.8e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
### blue = LARVAL & turquoise = LARVAL

# MALES:
res_M <- run.aov.mod(dat=datM,net=netM,phen=phenM)
kable(format(res_M,digits=2))
#|               |Larval  |Adult   |LarvalxAdult |fdr_Larval |fdr_Adult |fdr_LarvalxAdult |sign_Larval |sign_Adult |sign_LarvalxAdult |
#|:--------------|:-------|:-------|:------------|:----------|:---------|:----------------|:-----------|:----------|:-----------------|
#|MEturquoise    |6.5e-01 |3.2e-01 |2.2e-02      |7.2e-01    |4.4e-01   |7.3e-02          |0.0e+00     |0.0e+00    |0.0e+00           | 
#|MElightcyan    |3.0e-01 |4.1e-03 |7.0e-02      |4.4e-01    |2.1e-02   |1.6e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEgreen        |5.0e-10 |1.8e-03 |4.2e-02      |1.1e-08    |1.1e-02   |1.1e-01          |1.0e+00     |0.0e+00    |0.0e+00           | L
#|MEgreenyellow  |1.8e-03 |1.4e-06 |7.4e-01      |1.1e-02    |1.7e-05   |7.6e-01          |0.0e+00     |1.0e+00    |0.0e+00           | A
#|MElightyellow  |1.2e-02 |1.2e-01 |6.8e-01      |4.3e-02    |2.3e-01   |7.3e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEsalmon       |4.6e-03 |3.2e-01 |5.1e-01      |2.1e-02    |4.4e-01   |6.0e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEred          |3.3e-01 |9.3e-03 |6.3e-02      |4.4e-01    |3.6e-02   |1.5e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEyellow       |7.4e-02 |1.1e-01 |3.1e-01      |1.6e-01    |2.1e-01   |4.4e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEblack        |5.8e-02 |6.4e-01 |1.1e-01      |1.4e-01    |7.2e-01   |2.1e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEcyan         |1.8e-01 |1.7e-01 |7.8e-02      |3.0e-01    |2.9e-01   |1.6e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEroyalblue    |1.5e-02 |9.8e-13 |7.8e-01      |5.1e-02    |6.8e-11   |7.9e-01          |0.0e+00     |1.0e+00    |0.0e+00           | A
#|MEmidnightblue |6.7e-12 |1.6e-07 |6.6e-05      |2.3e-10    |2.8e-06   |6.5e-04          |1.0e+00     |1.0e+00    |1.0e+00           | L, A, LxA
#|MEpurple       |2.6e-02 |2.4e-02 |3.1e-01      |7.8e-02    |7.4e-02   |4.4e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEdarkred      |1.8e-03 |2.2e-01 |3.8e-01      |1.1e-02    |3.6e-01   |5.0e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MElightgreen   |5.1e-03 |2.6e-01 |8.1e-01      |2.2e-02    |4.1e-01   |8.1e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEdarkgreen    |2.7e-01 |8.0e-03 |7.2e-01      |4.1e-01    |3.3e-02   |7.6e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEgrey60       |6.5e-07 |3.4e-04 |3.3e-01      |8.9e-06    |2.9e-03   |4.4e-01          |1.0e+00     |0.0e+00    |0.0e+00           | L
#|MEpink         |4.5e-03 |1.4e-01 |1.3e-01      |2.1e-02    |2.4e-01   |2.4e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEbrown        |5.9e-04 |6.8e-01 |3.0e-02      |4.5e-03    |7.3e-01   |8.3e-02          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEblue         |4.2e-01 |2.7e-02 |4.3e-01      |5.3e-01    |7.9e-02   |5.3e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEmagenta      |4.5e-02 |6.1e-01 |4.9e-01      |1.2e-01    |7.1e-01   |6.0e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEtan          |5.2e-02 |2.5e-01 |6.0e-01      |1.3e-01    |4.0e-01   |7.0e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
#|MEgrey         |6.6e-01 |8.8e-02 |4.7e-01      |7.2e-01    |1.8e-01   |5.7e-01          |0.0e+00     |0.0e+00    |0.0e+00           |
### green = LARVAL, grey60 = LARVAL, greenyellow = ADULT, royalblue = ADULT, midnigthblue = LARVAL, ADULT, LARVAL x ADULT


