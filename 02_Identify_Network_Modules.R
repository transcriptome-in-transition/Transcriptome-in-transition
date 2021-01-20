####################################################################
### Identify network modules                                     ###
####################################################################

### Setting folders:
folderIN <- getwd()

folderIN_data <- paste0(folderIN,"/Workspaces")
pathIN_data <- paste0(folderIN_data,"/Normalized_data.Rdata")
folderIN_functions <- paste0(folderIN,"/RFunctions")
pathIN_functions <- paste0(folderIN_functions,c("/module.functions_V1.10.R","/DAVID.functions_V2.3.R"))
folderOUT_QC <- paste0(folderIN,"/QC")
tag <- "Day1_FB_17_07_2015"                  

### Loading packages & own util functions:
library("WGCNA")
library("cluster")
library("affy")

Dummi <- sapply(pathIN_functions,source)

### Loading data prepared data:
load(pathIN_data)

# FEMALES:
datF  <- exprs(Exgirl)
phenF <- pData(Exgirl)

## Remove outlier samples identified in first round of module making
## These were previously marked to spawn a module on their own ...
INDF  <- which(!(rownames(phenF) %in% c("7032-07-118.CEL","7032-07-282.CEL")))           
datF  <- datF[,INDF]
phenF <- phenF[INDF,]

# MALES:
datM  <- exprs(Exboy)
phenM <- pData(Exboy)
INDM  <- which(!(rownames(phenM) %in% c("7032-07-194.CEL")))                             
datM  <- datM[,INDM]
phenM <- phenM[INDM,]

### Check for whole-array outliers:
dir.create(folderOUT_QC)
pdf(paste0(folderOUT_QC,"/outliers_FEMALES_",tag,".pdf"),useDingbats=FALSE)
rec.eliminate.sample(datF,Zcut=2.5,main=NULL,plot=TRUE)
dev.off()
#=== Recursively eliminating outlier samples ===
#| Stopping elimination process at [0] with: 7032-07-236.CEL (Zk=-1.99)
#| Stopping elimination process at [1] with: 7032-07-236.CEL (Zk=-1.99)

pdf(paste0(folderOUT_QC,"/outliers_MALES_",tag,".pdf"),useDingbats=FALSE)
rec.eliminate.sample(datM,Zcut=2.5,main=NULL,plot=TRUE)
dev.off()
### 'VANILLA'
#| Eliminating [0]: 7032-07-206.CEL (Zk=-3.97)
#| Eliminating [1]: 7032-07-203.CEL (Zk=-3.02)
#| Eliminating [2]: 7032-07-114.CEL (Zk=-3.09)
#| Stopping elimination process at [3] with: 7032-07-100.CEL (Zk=-2.33) | 'no86': this becomes an outlier -2.64

## Remove outliers:
# FEMALES:
# MALES:
INDM <- which(!(colnames(datM) %in% c("7032-07-206.CEL","7032-07-203.CEL","7032-07-114.CEL")))                                      # 'VANILLA'
datM <- datM[,INDM]
phenM <- phenM[INDM,]
all(rownames(phenM)==colnames(datM))

### Plot sample dendrograms:
pdf(paste0(folderOUT_QC,"/sample_dendro_FEMALES_",tag,".pdf"),useDingbats=FALSE)
plot.sample.dendro(dat=datF,Zcut=2.5,pheno=phenF,main="Sample dendrogram and trait heatmap in FEMALES")
dev.off()
pdf(paste0(folderOUT_QC,"/sample_dendro_MALES_",tag,".pdf"),useDingbats=FALSE)
plot.sample.dendro(dat=datM,Zcut=2.5,pheno=phenM,main="Sample dendrogram and trait heatmap in MALES")
dev.off()

### Infer modules:
pathsOUT_net <- paste0(folderIN_data,"/netw_",c("FEMALES","MALES"),"_",tag,".RData")
netF <- comp.modules(datF,type="Tina")
save(netF,file=pathsOUT_net[1])
netM <- comp.modules(datM,type="Tina")
save(netM,file=pathsOUT_net[2])

### Plot module clustering:
pdf(paste0(folderOUT_QC,"/netw_FEMALES_",tag,".pdf"),useDingbats=FALSE,width=12)
plot.module.dendro(netF)
dev.off()
pdf(paste0(folderOUT_QC,"/netw_MALES_",tag,".pdf"),useDingbats=FALSE,width=12)
plot.module.dendro(netM)
dev.off()

### Check module-size:
pdf(paste0(folderOUT_QC,"/netw_FEMALES_sizes_",tag,".pdf"),useDingbats=FALSE,width=12)
plot.module.size(netF,main="Module Sizes in FEMALES")
dev.off()
pdf(paste0(folderOUT_QC,"/netw_MALES_sizes_",tag,".pdf"),useDingbats=FALSE,width=12)
plot.module.size(netM,main="Module Sizes in MALES")
dev.off()

### Plot module eigengene dendrograms:
pdf(paste0(folderOUT_QC,"/MEdendro_FEMALES_",tag,".pdf"),useDingbats=FALSE,width=12)
plot.ME.dendro(datF,netF)																# This does not work for only 3 modules ...
dev.off()
pdf(paste0(folderOUT_QC,"/MEdendro_MALES_",tag,".pdf"),useDingbats=FALSE,width=12)
plot.ME.dendro(datM,netM)
dev.off()

### Plot module eigengenes vs. sample similarities:
colF <- rep(rgb(0,0,1,0.5),nrow(phenF))
colF[which(phenM$Larval.food.level=="Control")] <- rgb(0,1,0,0.5)
colF[which(phenM$Larval.food.level=="Rich")] <- rgb(1,0,0,0.5)
pchF <- rep(22,nrow(phenM))
pchF[which(phenM$Adult.food.level=="C")] <- 21
pchF[which(phenM$Adult.food.level=="R")] <- 24
pdf(paste0(folderOUT_QC,"/MEvsZk_FEMALES_",tag,".pdf"),useDingbats=FALSE,width=12)
plot.sample.similarity.vs.module.eg(datF,netF,Zcut=2.5,hold=FALSE,colF,pchF)
dev.off()


colM <- rep(rgb(0,0,1,0.5),nrow(phenM))
colM[which(phenM$Larval.food.level=="Control")] <- rgb(0,1,0,0.5)
colM[which(phenM$Larval.food.level=="Rich")] <- rgb(1,0,0,0.5)
pchM <- rep(22,nrow(phenM))
pchM[which(phenM$Adult.food.level=="C")] <- 21
pchM[which(phenM$Adult.food.level=="R")] <- 24
pdf(paste0(folderOUT_QC,"/MEvsZk_MALES_",tag,".pdf"),useDingbats=FALSE,width=12)
plot.sample.similarity.vs.module.eg(datM,netM,Zcut=2.5,hold=FALSE,colM,pchM)
dev.off()
