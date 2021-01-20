#######################################################
### DAVID analyses                                  ###
#######################################################
library("RDAVIDWebService")

### Setting folders:
folderIN <- "."
folderIN_data <- paste0(folderIN,"/Workspaces")
pathIN_data <- paste0(folderIN_data,"/Normalized_data.Rdata")
folderIN_functions <- paste0(folderIN,"/RFunctions")
pathIN_functions <- paste0(folderIN_functions,c("/module.functions_V1.11.R","/DAVID.functions_V2.3.R"))
folderOUT_DAVID <- paste0(folderIN,"/David_output/Day1")
foldersOUT_DAVID <- paste0(folderOUT_DAVID,"/",c("FEMALES","MALES"))
tag <- "Day1"

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
pathsOUT_net <- paste0(folderIN_data,"/netw_",c("FEMALES","MALES"),"_",tag,".RData")
netF <- get(load("./Workspaces/netw_FEMALES_Day1_FB_17_07_2015.RData"))
netM <- get(load("./Workspaces/netw_MALES_Day1_FB_17_07_2015.RData"))

########################
### DAVID on MODULES ###
########################

david<-DAVIDWebService(email="chmay@bcgsc.ca", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
david
pathsOUT_DAVID <- paste0(folderOUT_DAVID,c("/DAVID_FEMALES_","/DAVID_MALES_"),tag,".RData")
names(pathsOUT_DAVID) <- c("FEMALES","MALES")
# FEMALES:
DAVID_F <- DAVID(net=netF,dat=datF,annot=c("GOTERM_BP_FAT","GOTERM_MF_FAT","GOTERM_CC_FAT"))
save(DAVID_F,file=pathsOUT_DAVID["FEMALES"])
# MALES:
DAVID_M <- DAVID(net=netM,dat=datM,annot=c("GOTERM_BP_FAT","GOTERM_MF_FAT","GOTERM_CC_FAT"))
save(DAVID_M,file=pathsOUT_DAVID["MALES"])
