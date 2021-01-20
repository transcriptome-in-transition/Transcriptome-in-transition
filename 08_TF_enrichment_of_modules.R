############################################################################
### 02_09_2015_Rerun_analysis_from_beginning: TFs Enrichments on Modules ###
############################################################################
library("WGCNA")
library("cluster")
library("affy")

### Setting folders:
folderIN <- "."
folderIN_data <- paste0(folderIN,"/Workspaces")
pathIN_data <- paste0(folderIN_data,"/Normalized_data.Rdata")
folderIN_functions <- paste0(folderIN,"/RFunctions")
pathIN_functions <- paste0(folderIN_functions,c("/module.functions_V1.11.R","/DAVID.functions_V2.3.R",
                                                "/PPI.functions_V1.R","/droid.functions_V1.2.R"))
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

#################################
### TF Enrichments on MODULES ###
#################################

# 1] Download a new version of TF-GENE mappings from DroID and store it locally (dropbox):
#dir.create("./DroID_results")
folderOUT_DroID <- paste0(folderIN,"/DroID_results")
#renew.DroID.TF(folderOUT_DroID)          # Preferably do this once per project ...
#=== Renewing DroID TF-Gene ===
# | Downloading data ...OK!
# | Mapping to Affy array ... 
#Removing  434  probes mapping to multiple genes 
#OK!
#| DroID TF-GENE mappings were stored under: 'DroID__TF__V2014_10__2015-09-02.RData'
#DONE!

# 2] Compute Enrichments:
res <- mod.droid.enr(dat=datM,net=netM,folderIN=folderOUT_DroID,type="TF")
#=== Performing TF Enrichments [DroID] === 
#| Loading prepared DroID TF-GENE Mappings: DroID__TF__V2014_10__2015-09-02.RData ... OK! 
#Removing  24  probes mapping to multiple genes 
#Removing  34  probes mapping to multiple genes 
#Removing  30  probes mapping to multiple genes 
#Removing  7  probes mapping to multiple genes 
#Removing  49  probes mapping to multiple genes 
#Removing  8  probes mapping to multiple genes 
#Removing  70  probes mapping to multiple genes 
#Removing  6  probes mapping to multiple genes 
#Removing  1  probes mapping to multiple genes 
#Removing  3  probes mapping to multiple genes 
#Removing  1  probes mapping to multiple genes 
#Removing  9  probes mapping to multiple genes 
#Removing  6  probes mapping to multiple genes 
#Removing  12  probes mapping to multiple genes 
#Removing  29  probes mapping to multiple genes 
#Removing  37  probes mapping to multiple genes 
#Removing  3  probes mapping to multiple genes 
#Removing  15  probes mapping to multiple genes 
#Removing  4  probes mapping to multiple genes 
#Removing  43  probes mapping to multiple genes 
#Removing  43  probes mapping to multiple genes
#etc.

# List topresults (enrichments only:)
# 1] Pruning results on numbers of annotations per TF:
TF <- names(which(res$Ngroup_test>=5))
# 2] Bonferroni correct:
minP <- 0.05/(ncol(res$p[,TF])*nrow(res$p[,TF]))
# 3] Select:
IND <- which((res$OR[,TF]>=1) & (res$p[,TF]<=minP),arr.ind=TRUE)
top <- data.frame(module=rownames(IND),TF=colnames(res$OR[,TF])[IND[,2]],OR=res$OR[,TF][IND],p=res$p[,TF][IND],
                  overlap=res$x11[,TF][IND],Ntf=res$x12[,TF][IND],Nmod=res$x21[,TF][IND],stringsAsFactors=FALSE)
knitr::kable(top[sort(top$p,index.return=TRUE,decreasing=FALSE)$ix,])
# Lots of blue ...

# Volcano plot of TF-GENE enrichments:
pdf(paste0(folderOUT_DroID,"/TF_MALES_volcano_",tag,".pdf"),useDingbats=FALSE,width=12)
plot.enr.volcano(res,groupnames=TF,main="Volcano plot of TF Enrichments")
dev.off()

# Heatmap:
STAT <- log2(res$OR[,TF])
par(mar = c(8, 8.5, 3, 3))
pdf(paste0(folderOUT_DroID,"/TF_MALES_heatmap_",tag,".pdf"),useDingbats=FALSE,width=12)
labeledHeatmap(Matrix=STAT,xLabels=colnames(STAT),yLabels=paste0("ME",rownames(STAT)),ySymbols=rownames(STAT),
               main="Modules vs TF",colors=blueWhiteRed(50),zlim=c(-5,5),cex.lab=0.8,font.lab=2)
dev.off()

####################################
### miRNA Enrichments on MODULES ###
####################################

# 1] Download a new version of miRNA-GENE mappings from DroID and store it locally (dropbox):
#renew.DroID.mir(folderOUT_DroID)           # Preferably do this once per project ...
#|=== Renewing DroID miRNA-Gene ===
#| Downloading data ...OK!
#| Mapping to Affy array ... 
#Removing  434  probes mapping to multiple genes 
#OK!
#| DroID miRNA-GENE mappings were stored under: 'DroID__mir__V2014_10__2015-09-03.RData'
#DONE!

# 2] Compute Enrichments:
res <- mod.droid.enr(dat=datM,net=netM,folderIN=folderOUT_DroID,type="MIR")

# List topresults (enrichments only:)
# 1] Pruning results on numbers of annotations per mir:
GROUP <- names(which(res$Ngroup_test>=5))
# 2] Bonferroni correct:
minP <- 0.05/(ncol(res$p[,GROUP])*nrow(res$p[,GROUP]))
# 3] Select:
IND <- which((res$OR[,GROUP]>=1) & (res$p[,GROUP]<=minP),arr.ind=TRUE)
top <- data.frame(module=rownames(IND),MIR=colnames(res$OR[,GROUP])[IND[,2]],OR=res$OR[,GROUP][IND],p=res$p[,GROUP][IND],
                  overlap=res$x11[,GROUP][IND],Ngroup=res$x12[,GROUP][IND],Nmod=res$x21[,GROUP][IND],stringsAsFactors=FALSE)
knitr::kable(top[sort(top$p,index.return=TRUE,decreasing=FALSE)$ix,])
# Lots of turqoise ...

# Volcano plot of TF-GENE enrichments:
pdf(paste0(folderOUT_DroID,"/MIR_MALES_volcano_",tag,".pdf"),useDingbats=FALSE,width=12)
plot.enr.volcano(res,groupnames=GROUP,main="Volcano plot of MIR Enrichments")
dev.off()

# Heatmap:
STAT <- log2(res$OR[,GROUP])
par(mar = c(8, 8.5, 3, 3))
pdf(paste0(folderOUT_DroID,"/MIR_MALES_heatmap_",tag,".pdf"),useDingbats=FALSE,width=20)
labeledHeatmap(Matrix=STAT,xLabels=colnames(STAT),yLabels=paste0("ME",rownames(STAT)),ySymbols=rownames(STAT),
               main="Modules vs MIR",colors=blueWhiteRed(50),zlim=c(-5,5),cex.lab=0.8,font.lab=2)
dev.off()



