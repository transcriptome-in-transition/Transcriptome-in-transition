##Functions to plot eigengenes and to plot median expression for tissues/treatments etc in order of male dendrogram
#V.1.1: Reorder plots of FlyAtlas Median Expression by tissue clustering
##Plot by age
mod.tis.spec.age <- function(dat,net,pathIN_TA,plot=FALSE,median=TRUE){
  cat("=== Performing Tissue Specificity Analyses === \n")
  ## Load prepared Tissue Atlas data:
  cat(paste("| Preparing Tissue Atlas Data ... "))
  FA <- get(load(pathIN_TA))
  pheno_FA <- as(phenoData(FA),"data.frame")
  dat_FA <- exprs(FA)
  ## Compute medians per Module per Tissue in TA:
  if(median){
    mods <- sapply(get.mod.members(net,dat),names)
    if(length(mods) > 3) mods<-mods[c(22,11,7,8,13,20,18,23,1,4,19,15,17,6,12,5,10,16,3,2,14,21,9)]
    tissIND <- split(1:nrow(pheno_FA),droplevels(as.factor(pheno_FA$Time.point)))
    MED <- sapply(mods,function(x) sapply(tissIND,function(y) median(as.vector(dat_FA[x,y]),na.rm=TRUE)))
    if(plot){
      labeledHeatmap(t(MED),xLabels=names(tissIND),colors = blueWhiteRed(50),yLabels=paste0("ME",names(mods)),ySymbols=names(mods),
                     main="Median EigenGene Expression by Age",cex.lab.x=0.5,invertColors=TRUE)
    }
    return(MED)
  }  
}

###Plot by larval diet##
mod.tis.spec.lfl <- function(dat,net,pathIN_TA,plot=FALSE,median=TRUE){
  cat("=== Performing Tissue Specificity Analyses === \n")
  ## Load prepared Tissue Atlas data:
  cat(paste("| Preparing Tissue Atlas Data ... "))
  FA <- get(load(pathIN_TA))
  pheno_FA <- as(phenoData(FA),"data.frame")
  dat_FA <- exprs(FA)
  ## Compute medians per Module per Tissue in TA:
  if(median){
    mods <- sapply(get.mod.members(net,dat),names)
    if(length(mods) > 3) mods<-mods[c(22,11,7,8,13,20,18,23,1,4,19,15,17,6,12,5,10,16,3,2,14,21,9)]
    tissIND <- split(1:nrow(pheno_FA),droplevels(as.factor(pheno_FA$Larval.food.level)))
    MED <- sapply(mods,function(x) sapply(tissIND,function(y) median(as.vector(dat_FA[x,y]),na.rm=TRUE)))
    if(plot){
      labeledHeatmap(t(MED),xLabels=names(tissIND),colors = blueWhiteRed(50),yLabels=paste0("ME",names(mods)),ySymbols=names(mods),
                     main="Median EigenGene Expression by larval diet",cex.lab.x=0.5,invertColors=TRUE)
    }
    return(MED)
  }
}  

###Plot by adult diet##
mod.tis.spec.afl <- function(dat,net,pathIN_TA,plot=FALSE,median=TRUE){
  cat("=== Performing Tissue Specificity Analyses === \n")
  ## Load prepared Tissue Atlas data:
  cat(paste("| Preparing Tissue Atlas Data ... "))
  FA <- get(load(pathIN_TA))
  pheno_FA <- as(phenoData(FA),"data.frame")
  dat_FA <- exprs(FA)
  ## Compute medians per Module per Tissue in TA:
  if(median){
    mods <- sapply(get.mod.members(net,dat),names)
    if(length(mods) > 3) mods<-mods[c(22,11,7,8,13,20,18,23,1,4,19,15,17,6,12,5,10,16,3,2,14,21,9)]
    tissIND <- split(1:nrow(pheno_FA),droplevels(as.factor(pheno_FA$Adult.food.level)))
    MED <- sapply(mods,function(x) sapply(tissIND,function(y) median(as.vector(dat_FA[x,y]),na.rm=TRUE)))
    if(plot){
      labeledHeatmap(t(MED),xLabels=names(tissIND),colors = blueWhiteRed(50),yLabels=paste0("ME",names(mods)),ySymbols=names(mods),
                     main="Median EigenGene Expression by adult diet",cex.lab.x=0.5,invertColors=TRUE)
    }
    return(MED)
  }
}  

##Plot by age, larval diet, and adult diet
mod.tis.spec.TAL <- function(dat,net,pathIN_TA,plot=FALSE,median=TRUE){
  cat("=== Performing Tissue Specificity Analyses === \n")
  ## Load prepared Tissue Atlas data:
  cat(paste("| Preparing Tissue Atlas Data ... "))
  FA <- get(load(pathIN_TA))
  pheno_FA <- as(phenoData(FA),"data.frame")
  dat_FA <- exprs(FA)
  ## Compute medians per Module per Tissue in TA:
  if(median){
    mods <- sapply(get.mod.members(net,dat),names)
    if(length(mods) > 3) mods<-mods[c(22,11,7,8,13,20,18,23,1,4,19,15,17,6,12,5,10,16,3,2,14,21,9)]
    tissIND <- split(1:nrow(pheno_FA),droplevels(as.factor(pheno_FA$TAL)))
    MED <- sapply(mods,function(x) sapply(tissIND,function(y) median(as.vector(dat_FA[x,y]),na.rm=TRUE)))
    if(plot){
      labeledHeatmap(t(MED),xLabels=names(tissIND),colors = blueWhiteRed(50),yLabels=paste0("ME",names(mods)),ySymbols=names(mods),
                     main="Median EigenGene Expression by TAL",cex.lab.x=0.5,invertColors=TRUE)
    }
    return(MED)
  }
}

##Function to plot median expression of our modules in FlyAtlas Data##
mod.tis.spec.tina <- function(dat,net,pathIN_TA,plot=FALSE,median=TRUE){
  cat("=== Performing Tissue Specificity Analyses === \n")
  ## Load prepared Tissue Atlas data:
  cat(paste("| Preparing Tissue Atlas Data ... "))
  FA <- get(load(pathIN_TA))
  pheno_FA <- as(phenoData(FA),"data.frame")
  dat_FA <- exprs(FA)
  ## Compute medians per Module per Tissue in TA:
  if(median){
    mods <- sapply(get.mod.members(net,dat),names)
    if(length(mods) > 3) mods<-mods[c(22,11,7,8,13,20,18,23,1,4,19,15,17,6,12,5,10,16,3,2,14,21,9)]
    tissIND <- split(1:nrow(pheno_FA),droplevels(pheno_FA$tissue))
    tissIND<-tissIND[c("Adult_Testes","Adult_Ovary","Adult_Wings","Adult_Brain","Adult_Thoracoabdominal_ganglion","Adult_Eye","Adult_Trachea","Adult_Hind_Gut","Adult_Mid_Gut","Adult_Thoracic_Muscle","Adult_Accessory_gland","Adult_Salivary_Gland","Adult_Crop","Adult_Fatbody","Adult_Heart","Adult_Male_Ejaculatory_Duct")]
    MED <- sapply(mods,function(x) sapply(tissIND,function(y) median(as.vector(dat_FA[x,y]),na.rm=TRUE)))
    if(plot){
      labeledHeatmap(t(MED),xLabels=names(tissIND),colors = blueWhiteRed(50),yLabels=paste0("ME",names(mods)),ySymbols=names(mods),
                     main="Median Expression in Fly Atlas",cex.lab.x=0.5,invertColors=TRUE)
    }
    return(MED)
  }
}  

mod.tis.spec.tina.key <- function(dat,net,pathIN_TA,plot=FALSE,median=TRUE){
        cat("=== Performing Tissue Specificity Analyses === \n")
        ## Load prepared Tissue Atlas data:
        cat(paste("| Preparing Tissue Atlas Data ... "))
        FA <- get(load(pathIN_TA))
        pheno_FA <- as(phenoData(FA),"data.frame")
        dat_FA <- exprs(FA)
        ## Compute medians per Module per Tissue in TA:
        if(median){
                mods <- sapply(get.mod.members(net,dat),names)
                mods<-mods[c(7,8,19,15,10)]
                #if(length(mods) > 3) mods<-mods[c(22,11,7,8,13,20,18,23,1,4,19,15,17,6,12,5,10,16,3,2,14,21,9)]
                tissIND <- split(1:nrow(pheno_FA),droplevels(pheno_FA$tissue))
                tissIND<-tissIND[c("Adult_Testes","Adult_Ovary","Adult_Wings","Adult_Brain","Adult_Thoracoabdominal_ganglion","Adult_Eye","Adult_Trachea","Adult_Hind_Gut","Adult_Mid_Gut","Adult_Thoracic_Muscle","Adult_Accessory_gland","Adult_Salivary_Gland","Adult_Crop","Adult_Fatbody","Adult_Heart","Adult_Male_Ejaculatory_Duct")]
                MED <- sapply(mods,function(x) sapply(tissIND,function(y) median(as.vector(dat_FA[x,y]),na.rm=TRUE)))
                if(plot){
                        labeledHeatmap(t(MED),xLabels=names(tissIND),colors = blueWhiteRed(50),yLabels=paste0("ME",names(mods)),ySymbols=names(mods),
                                       main="Median Expression in Fly Atlas",cex.lab.x=0.5,invertColors=TRUE)
                }
                return(MED)
        }
}  


mod.tis.spec.gut <- function(dat,net,ExSet,plot=FALSE,median=TRUE){
  cat("=== Performing Tissue Specificity Analyses === \n")
  ## Load prepared Tissue Atlas data:
  cat(paste("| Preparing Tissue Atlas Data ... "))
  pheno_FA <- as(phenoData(ExSet),"data.frame")
  dat_FA <- exprs(ExSet)
  ## Compute medians per Module per Tissue in TA:
  if(median){
    mods <- sapply(get.mod.members(net,dat),names)
    if(length(mods) > 3) mods<-mods[c(22,11,7,8,13,20,18,23,1,4,19,15,17,6,12,5,10,16,3,2,14,21,9)]
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

mod.tis.spec.gut.only <- function(dat,net,ExSet,plot=FALSE,median=TRUE){
  cat("=== Performing Tissue Specificity Analyses === \n")
  ## Load prepared Tissue Atlas data:
  cat(paste("| Preparing Tissue Atlas Data ... "))
  pheno_FA <- as(phenoData(ExSet),"data.frame")
  dat_FA <- exprs(ExSet)
  ## Compute medians per Module per Tissue in TA:
  if(median){
    mods <- sapply(get.mod.members(net,dat),names)
    mods<-mods[c("royalblue", "midnightblue")]
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


mod.tis.spec.tina.larval <- function(dat,net,pathIN_TA,plot=FALSE,median=TRUE){
  cat("=== Performing Tissue Specificity Analyses === \n")
  ## Load prepared Tissue Atlas data:
  cat(paste("| Preparing Tissue Atlas Data ... "))
  ## Compute medians per Module per Tissue in TA:
  if(median){
    mods <- sapply(get.mod.members(net,dat),names)
    if(length(mods) > 3) mods<-mods[c(22,11,7,8,13,20,18,23,1,4,19,15,17,6,12,5,10,16,3,2,14,21,9)]
    tissIND <- split(1:nrow(pheno_FA),droplevels(pheno_FA$tissue))
    tissIND<-tissIND[c("Larval_Wandering_fat_body","Larval_Feeding_Central_Nevous_System","Larval_Feeding_Mid_Gut","Larvae_Wandering_Tubules","Larval_Feeding_Malpighian_Tubule","Larval_Feeding_Hind_Gut","Whole_Larvae_Feeding","Larval_Feeding_Carcass","Larval_Feeding_Trachea","Larval_Feeding_Fatbody","Larval_Feeding_Salivary_Gland")]
    MED <- sapply(mods,function(x) sapply(tissIND,function(y) median(as.vector(dat_FA[x,y]),na.rm=TRUE)))
    if(plot){
      labeledHeatmap(t(MED),xLabels=names(tissIND),colors = blueWhiteRed(50),yLabels=paste0("ME",names(mods)),ySymbols=names(mods),
                     main="Median Expression in Larval Tissues in Fly Atlas",cex.lab.x=0.5,invertColors=TRUE)
    }
    return(MED)
  }
}  

##############################################
###Plot at later time points only key modules
##Plot by age
mod.tis.spec.age.key <- function(dat,net,pathIN_TA,plot=FALSE,median=TRUE){
  cat("=== Performing Tissue Specificity Analyses === \n")
  ## Load prepared Tissue Atlas data:
  cat(paste("| Preparing Tissue Atlas Data ... "))
  FA <- get(load(pathIN_TA))
  pheno_FA <- as(phenoData(FA),"data.frame")
  dat_FA <- exprs(FA)
  ## Compute medians per Module per Tissue in TA:
  if(median){
    mods <- sapply(get.mod.members(net,dat),names)
    if(length(mods) > 3) mods<-mods[c("green","greenyellow","royalblue","midnightblue","grey60")]
    tissIND <- split(1:nrow(pheno_FA),droplevels(as.factor(pheno_FA$Time.point)))
    MED <- sapply(mods,function(x) sapply(tissIND,function(y) median(as.vector(dat_FA[x,y]),na.rm=TRUE)))
    if(plot){
      labeledHeatmap(t(MED),xLabels=names(tissIND),colors = blueWhiteRed(50),yLabels=paste0("ME",names(mods)),ySymbols=names(mods),
                     main="Median EigenGene Expression by Age",cex.lab.x=0.5,invertColors=TRUE)
    }
    return(MED)
  }  
}

###Plot by larval diet##
mod.tis.spec.lfl.key <- function(dat,net,pathIN_TA,plot=FALSE,median=TRUE){
  cat("=== Performing Tissue Specificity Analyses === \n")
  ## Load prepared Tissue Atlas data:
  cat(paste("| Preparing Tissue Atlas Data ... "))
  FA <- get(load(pathIN_TA))
  pheno_FA <- as(phenoData(FA),"data.frame")
  dat_FA <- exprs(FA)
  ## Compute medians per Module per Tissue in TA:
  if(median){
    mods <- sapply(get.mod.members(net,dat),names)
    if(length(mods) > 3) mods<-mods[c("green","greenyellow","royalblue","midnightblue","grey60")]
    tissIND <- split(1:nrow(pheno_FA),droplevels(as.factor(pheno_FA$Larval.food.level)))
    MED <- sapply(mods,function(x) sapply(tissIND,function(y) median(as.vector(dat_FA[x,y]),na.rm=TRUE)))
    if(plot){
      labeledHeatmap(t(MED),xLabels=names(tissIND),colors = blueWhiteRed(50),yLabels=paste0("ME",names(mods)),ySymbols=names(mods),
                     main="Median EigenGene Expression by larval diet",cex.lab.x=0.5,invertColors=TRUE)
    }
    return(MED)
  }
}  

###Plot by adult diet##
mod.tis.spec.afl.key <- function(dat,net,pathIN_TA,plot=FALSE,median=TRUE){
  cat("=== Performing Tissue Specificity Analyses === \n")
  ## Load prepared Tissue Atlas data:
  cat(paste("| Preparing Tissue Atlas Data ... "))
  FA <- get(load(pathIN_TA))
  pheno_FA <- as(phenoData(FA),"data.frame")
  dat_FA <- exprs(FA)
  ## Compute medians per Module per Tissue in TA:
  if(median){
    mods <- sapply(get.mod.members(net,dat),names)
    if(length(mods) > 3) mods<-mods[c("green","greenyellow","royalblue","midnightblue","grey60")]
    tissIND <- split(1:nrow(pheno_FA),droplevels(as.factor(pheno_FA$Adult.food.level)))
    MED <- sapply(mods,function(x) sapply(tissIND,function(y) median(as.vector(dat_FA[x,y]),na.rm=TRUE)))
    if(plot){
      labeledHeatmap(t(MED),xLabels=names(tissIND),colors = blueWhiteRed(50),yLabels=paste0("ME",names(mods)),ySymbols=names(mods),
                     main="Median EigenGene Expression by adult diet",cex.lab.x=0.5,invertColors=TRUE)
    }
    return(MED)
  }
}  

##Plot by age, larval diet, and adult diet
mod.tis.spec.TAL.key <- function(dat,net,pathIN_TA,plot=FALSE,median=TRUE){
  cat("=== Performing Tissue Specificity Analyses === \n")
  ## Load prepared Tissue Atlas data:
  cat(paste("| Preparing Tissue Atlas Data ... "))
  FA <- get(load(pathIN_TA))
  pheno_FA <- as(phenoData(FA),"data.frame")
  dat_FA <- exprs(FA)
  ## Compute medians per Module per Tissue in TA:
  if(median){
    mods <- sapply(get.mod.members(net,dat),names)
    if(length(mods) > 3) mods<-mods[c("green","greenyellow","royalblue","midnightblue","grey60")]
    tissIND <- split(1:nrow(pheno_FA),droplevels(as.factor(pheno_FA$TAL)))
    MED <- sapply(mods,function(x) sapply(tissIND,function(y) median(as.vector(dat_FA[x,y]),na.rm=TRUE)))
    if(plot){
      labeledHeatmap(t(MED),xLabels=names(tissIND),colors = blueWhiteRed(50),yLabels=paste0("ME",names(mods)),ySymbols=names(mods),
                     main="Median EigenGene Expression by TAL",cex.lab.x=0.5,invertColors=TRUE)
    }
    return(MED)
  }
}



