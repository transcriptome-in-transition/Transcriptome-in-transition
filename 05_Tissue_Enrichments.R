####################################################################
###                 Tissue Enrichments                           ###
####################################################################
library("WGCNA")
library("cluster")
library("affy")

### Setting folders:
folderIN <- "."
folderIN_data <- paste0(folderIN,"/Workspaces")
pathIN_data <- paste0(folderIN_data,"/Normalized_data.Rdata")
folderIN_functions <- paste0(folderIN,"/RFunctions")
pathIN_functions <- paste0(folderIN_functions,c("/module.functions_V1.8.R","/DAVID.functions_V2.3.R", "/eigengene_and_median_functions_V1.0.R"))
folderOUT_QC <- paste0(folderIN,"/QC")
folderOUT_TISSUE <- paste0(folderIN,"/Tissue_output")
tag <- "Day1_FB_17_07_2015"

### Load data:
load(pathIN_data)
# Phenotypes:
phenM <- phenoData(Exboy)
phenF <- phenoData(Exgirl)
# Data:
datM <- exprs(Exboy)
datF <- exprs(Exgirl)

### Loading inferred networks:
load("./Workspaces/netw_FEMALES_Day1_FB_17_07_2015.RData")
load("./Workspaces/netw_MALES_Day1_FB_17_07_2015.RData")

### Loading own util functions:
Dummi <- sapply(pathIN_functions,source)

###############################################
### Do some diagnostics on inferred modules ###
###############################################

# 1] Check module-sizes:
pdf(paste0(folderOUT_QC,"/netw_MALES_sizes_",tag,".pdf"),useDingbats=FALSE,width=12,height=12)
plot.module.size(netM,main="Module Sizes in MALES")
dev.off()
pdf(paste0(folderOUT_QC,"/netw_FEMALES_sizes_",tag,".pdf"),useDingbats=FALSE,width=12,height=12)
plot.module.size(netF,main="Module Sizes in FEMALES")
dev.off()

###################################################
### Tissue Enrichments on Expression Abundances ###
###################################################
pathIN_FA <- "./External_datasets/FlyAtlas/Expr_adult_FA_07_05_2015.RData"
enrM_exp <- mod.tis.enr(netM,datM,pathIN_TA=pathIN_FA,Ecut=5,background="expressed")
#=== Performing Tissue Enrichments === 
#| Preparing Tissue Atlas Data ... OK! 
#| Setting background : 'expressed' [N = 12400] 
#| Testing for enrichments per tissue at Absolute Expression >= 5...| Selecting genes expressed in ANY tissue of the atlas ...
#|   Adult_Accessory_gland [# Expr >= 5: 6741]: 
#******************************
#|   Adult_Brain [# Expr >= 5: 7570]: 
#******************************
#|   Adult_Crop [# Expr >= 5: 7152]: 
#******************************
#|   Adult_Eye [# Expr >= 5: 7572]: 
#******************************
#|   Adult_Fatbody [# Expr >= 5: 6778]: 
#******************************
#|   Adult_Heart [# Expr >= 5: 7191]: 
#******************************
#|   Adult_Hind_Gut [# Expr >= 5: 7389]: 
#******************************
#|   Adult_Male_Ejaculatory_Duct [# Expr >= 5: 6965]: 
#******************************
#|   Adult_Mid_Gut [# Expr >= 5: 6967]: 
#******************************
#|   Adult_Ovary [# Expr >= 5: 6998]: 
#******************************
#|   Adult_Salivary_Gland [# Expr >= 5: 6407]: 
#******************************
#|   Adult_Testes [# Expr >= 5: 7696]: 
#******************************
#|   Adult_Thoracic_Muscle [# Expr >= 5: 6418]: 
#******************************
#|   Adult_Thoracoabdominal_ganglion [# Expr >= 5: 7457]: 
#******************************
#|   Adult_Trachea [# Expr >= 5: 7487]: 
#******************************
#|   Adult_Wings [# Expr >= 5: 7489]: 
#******************************
# === DONE! ===
dir.create("./Tissue_output")
pdf(paste0(folderOUT_TISSUE,"/TISSUE_MALES_expr_adult_",tag,".pdf"),useDingbats=FALSE,width=12)
plot.tis.enr(enrM_exp)
dev.off()

mod.tis.spec.tina(datM,netM,pathIN_FA,plot=TRUE,median=TRUE)


enrF_exp <- mod.tis.enr(netF,datF,pathIN_TA=pathIN_FA,Ecut=5,background="expressed")
#=== Performing Tissue Enrichments === 
#| Preparing Tissue Atlas Data ... OK! 
#| Setting background : 'expressed' [N = 12400] 
#| Testing for enrichments per tissue at Absolute Expression >= 5...| Selecting genes expressed in ANY tissue of the atlas ...
#|   Adult_Accessory_gland [# Expr >= 5: 6741]: 
#******
#|   Adult_Brain [# Expr >= 5: 7570]: 
#******
#|   Adult_Crop [# Expr >= 5: 7152]: 
#******
#|   Adult_Eye [# Expr >= 5: 7572]: 
#******
#|   Adult_Fatbody [# Expr >= 5: 6778]: 
#******
#|   Adult_Heart [# Expr >= 5: 7191]: 
#******
#|   Adult_Hind_Gut [# Expr >= 5: 7389]: 
#******
#|   Adult_Male_Ejaculatory_Duct [# Expr >= 5: 6965]: 
#******
#|   Adult_Mid_Gut [# Expr >= 5: 6967]: 
#******
#|   Adult_Ovary [# Expr >= 5: 6998]: 
#******
#|   Adult_Salivary_Gland [# Expr >= 5: 6407]: 
#******
#|   Adult_Testes [# Expr >= 5: 7696]: 
#******
#|   Adult_Thoracic_Muscle [# Expr >= 5: 6418]: 
#******
#|   Adult_Thoracoabdominal_ganglion [# Expr >= 5: 7457]: 
#******
#|   Adult_Trachea [# Expr >= 5: 7487]: 
#******
#|   Adult_Wings [# Expr >= 5: 7489]: 
#******
# === DONE! ===
pdf(paste0(folderOUT_TISSUE,"/TISSUE_FEMALES_expr_adult_",tag,".pdf"),useDingbats=FALSE,width=12)
plot.tis.enr(enrF_exp)
dev.off()

mod.tis.spec.tina(datF,netF,pathIN_FA,plot=TRUE,median=TRUE)

### Plot expression in gut
##Read in gut dataset
files<-("./Gut_region_specificity")
setwd(files)

##Read in CELfiles##
CELfiles<-list.celfiles()
CELfiles
pd <- read.AnnotatedDataFrame("./sample_data.txt",header=TRUE,row.names=1)
rawdata <- ReadAffy(filenames=CELfiles,phenoData=pd,verbose=TRUE)

GutExprs<-rma(rawdata)
pdf("Median Expression in Gut Tissues of all modules (Males).pdf",useDingbats=FALSE,width=12)
mod.tis.spec.gut(datM,netM,GutExprs,plot=TRUE,median=TRUE)
dev.off()

pdf("Median Expression in Gut Tissues of all modules (Females).pdf",useDingbats=FALSE,width=12)
mod.tis.spec.gut(datF,netF,GutExprs,plot=TRUE,median=TRUE)
dev.off()

pdf("Median Expression in Gut Tissues from Fly Gut Database (Males).pdf",useDingbats=FALSE,width=12)
mod.tis.spec.gut.only(datM,netM,GutExprs,plot=TRUE,median=TRUE)
dev.off()

mod.tis.spec.gut.top.five <- function(dat,net,ExSet,plot=FALSE,median=TRUE){
  cat("=== Performing Tissue Specificity Analyses === \n")
  ## Load prepared Tissue Atlas data:
  cat(paste("| Preparing Tissue Atlas Data ... "))
  pheno_FA <- as(phenoData(ExSet),"data.frame")
  dat_FA <- exprs(ExSet)
  ## Compute medians per Module per Tissue in TA:
  if(median){
    mods <- sapply(get.mod.members(net,dat),names)
    if(length(mods) > 3) if(length(mods) > 3) mods<-mods[c("green","greenyellow","royalblue","midnightblue","grey60")]
    tissIND <- split(1:nrow(pheno_FA),droplevels(as.factor(pheno_FA$Characteristics.phenotype.)))
    tissIND<-tissIND[c("Crop","proventriculus/R1","R2","R3","R4","R5","Hindgut")]
    MED <- sapply(mods,function(x) sapply(tissIND,function(y) median(as.vector(dat_FA[x,y]),na.rm=TRUE)))
    if(plot){
      labeledHeatmap(t(MED),xLabels=names(tissIND),colors = blueWhiteRed(50),yLabels=paste0("ME",names(mods)),ySymbols=names(mods),
                     main="Median Expression in Gut Tissues from Fly Gut Database",cex.lab.x=0.5,invertColors=TRUE)
    }
    return(MED)
  }
}  
pdf("Median Expression in Gut Tissues of top five modules (Males).pdf",useDingbats=FALSE,width=12)
mod.tis.spec.gut.top.five(datM,netM,GutExprs,plot=TRUE,median=TRUE)
dev.off()



