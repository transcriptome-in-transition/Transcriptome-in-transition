
## DroID.functions:

# Aim: To handle DroID data ...

# Change-log:
# V1.0	03-09-2015	Erik	Initiation of this file

# REMARKS:
# Let's store local copies of DroID resources in 
# '~/15_07_2015_Rerun_analysis_from_beginning/DroID_results' under:
# DroID__dataID__tag__timestamp.RData
# Let's then define functions:
# 'get.DroID.dataID':   To download the resource
# 'renew.DroID.dataID': To make a new fresh local time-stamped copy
# 'load.DroID.dataID':  To load the latest local copy/dowload the latest version in your workspace 

## Defining a SubFunction to timestamp my dowloads:
timestamp <- function(digits=0,dateOnly=TRUE){
  old <- getOption("digits.secs")
  options(digits.secs = digits)
  stamp <- gsub("[. :]","_",Sys.time())
  if(dateOnly){
    stamp <- unlist(strsplit(stamp,"_"))[1]
  }
  options(digits.secs = old)
  return(stamp)
}

## Defining a SubFunction to construct consistent paths:
make.path <- function(folderOUT=NULL,dataID,tag,db="DroID"){
  if(!is.null(folderOUT)){
    pathOUT <- paste0(folderOUT,"/",db,"__",dataID,"__",tag,"__",timestamp(),".RData")
  } else {
    pathOUT <- paste0(db,"__",dataID,"__",tag,"__",timestamp(),".RData")
  }
  return(pathOUT)
}

## Defining a SubFunction for returing the path to the latest version:
get.path <- function(folderIN,dataID,db=NULL,latest=TRUE){
  paths <- dir(folderIN,full.names=TRUE)
  # 1) Subselect on dataID:
  paths <- paths[grep(paste0("__",dataID,"__"),basename(paths))]
  # 2) Subselect on database (if defined):
  if(!is.null(db)){
    paths <- paths[grep(paste0(db,"__"),basename(paths))]
  }
  # 3) Extract dates to find latest entry:
  if(latest & length(paths)>1){
    dates <- sapply(strsplit(basename(paths),"[.]|__"),function(x) grep("2[0-9]{3}-[0-9]{2}-[0-9]{2}",unlist(x),value=TRUE))
    paths <- paths[sort(dates,decreasing=TRUE,index.return=TRUE)$ix[1]]
  }
  return(paths)
}

### TF-GENE:

## Defining a SubFunction for downloading TF-GENE links from DroID:
get.DroID.TF <- function(pathIN="http://www.droidb.org/data/DroID_v2014_10/tf_gene.txt"){
  # 1] Downloading data:
  cat("| Downloading data ...")
  dat <- read.delim(pathIN,stringsAsFactors = FALSE)
  cat("OK!\n")
  # 2] Append AffyIDs:
  cat("| Mapping to Affy array ... \n")
  FB2Affy <- Affy2FlyBase(affyIDs=NULL,rev=TRUE)
  FB2Affy <- data.frame(fbID=names(FB2Affy),affyID=as.vector(FB2Affy),stringsAsFactors = FALSE)
  dat <- merge(dat,FB2Affy,by.x="FLY_TARGET_GENE",by.y="fbID",all.x=TRUE)
  colnames(dat)[which(colnames(dat)=="affyID")] <- "affyID_GENE"
  dat <- merge(dat,FB2Affy,by.x="FLY_TF_GENE",by.y="fbID",all.x=TRUE)
  colnames(dat)[which(colnames(dat)=="affyID")] <- "affyID_TF"
  cat("OK!\n")
  return(dat)
}

## Defining a SubFunction for storing a fresh copy of TF-Gene links from DroID:
renew.DroID.TF <- function(folderOUT,pathIN="http://www.droidb.org/data/DroID_v2014_10/tf_gene.txt"){
  cat("=== Renewing DroID TF-Gene ===\n")
  # 1) Download & map to affyIDs:
  TF <- get.DroID.TF(pathIN=pathIN)
  # 2) Store:
  tag <- toupper(gsub("DroID_","",grep("DroID",unlist(strsplit(pathIN,"/")),value=TRUE)))
  pathOUT <- make.path(folderOUT=folderOUT,dataID="TF",tag=tag)
  save(TF,file=pathOUT)
  cat(paste0("| DroID TF-GENE mappings were stored under: '",basename(pathOUT),"'\n"))
  cat("DONE!\n")
}

## Defining a SubFunction for loading the latest version of TF-GENE:
load.DroID.TF <- function(folderIN){
  pathIN <- get.path(folderIN,dataID="TF",db="DroID",latest=TRUE)
  if(length(pathIN)==0){
    stop("Could not find a local copy in the designated path ...")
  }
  cat(paste0("| Loading prepared DroID TF-GENE Mappings: ",basename(pathIN)," ..."))
  TF <- get(load(pathIN))
  cat(" OK! \n")
  return(TF)
}

## Defining a SubFunction for performing a 2x2 adjusted Fisher exact test:
buildandtest2by2 <- function(X,Y,ALL){
  cat("*")
  if(is.null(X)){
    X <- ""
  }
  x11 <- length(intersect(X,Y))
  x12 <- length(setdiff(X,x11))
  x21 <- length(setdiff(Y,x11))
  x22 <- length(setdiff(ALL,x11))
  # Test:
  res <- fisher.test(matrix(c(x11,x12,x21,x22),nrow=2,byrow=TRUE))
  return(list(p=res$p.value,OR=res$estimate[[1]],low=res$conf.int[1],high=res$conf.int[2],x11=x11,x12=x12,x21=x21,x22=x22))
}

## Defining a SubFunction for performing an ontology-prune prior to 2x2 testing:
# Similar as to what DAVID does, you can prune your foreground and background to all genes known to a 
# certain ontology ('system'). Example: If testing for the overlap of your genelist with lists of genes 
# known to be expressed in a tissue, you may want to consider genes that are expressed at least once in
# ANY tissue, especially if the number of tissues is large. If not, you bias the background towards 
# 'rubbish', hence, you are likely to find enrichments.
ontologyprune <- function(geneset1,geneset2){
  # Prune doubles and missing annotations from genesets:
  geneset1 <- sapply(geneset1,function(x) setdiff(unique(na.omit(x)),""))
  geneset2 <- sapply(geneset2,function(x) setdiff(unique(na.omit(x)),""))
  # Count:
  N1_1 <- sapply(geneset1,length)
  N1tot_1 <- length(unique(unlist(geneset1)))
  N1_2 <- sapply(geneset2,length)
  N1tot_2 <- length(unique(unlist(geneset2)))
  # Prune:
  background <- intersect(unlist(geneset1),unlist(geneset2))
  geneset1 <- sapply(geneset1,function(x) intersect(x,background))
  geneset2 <- sapply(geneset2,function(x) intersect(x,background))
  # Count:
  N2_1 <- sapply(geneset1,length)
  N2tot_1 <- length(unique(unlist(geneset1)))
  N2_2 <- sapply(geneset2,length)
  N2tot_2 <- length(unique(unlist(geneset2)))
  return(list(geneset1=geneset1,geneset2=geneset2,background=background,N1_start=N1_1,N1tot_start=N1tot_1,
              N2_start=N1_2,N2tot_start=N1tot_2,N1_end=N2_1,N1tot_end=N2tot_1,N2_end=N2_2,N2tot_end=N2tot_2))
}

## Defining a SubFunction for computing a TF enrichment for modules:
mod.tf_droid.enr <- function(dat,net,folderIN=NULL){
  cat("=== Performing TF Enrichments [DroID] === \n")
  if(is.null(folderIN)){
    TF <- get.DroID.TF()
  } else {
    TF <- load.DroID.TF(folderIN)
  }
  TF <- unique(TF[,c("FLY_TARGET_GENE","FLY_TF_GENE","GENE_SYMBOL","TF_SYMBOL","affyID_GENE")])
  # Extracting modules and map them to genes:
  geneList  <- sapply(sapply(get.mod.members(net,dat),names),function(x) setdiff(unique(na.omit(Affy2FlyBase(x))),""))
  # Extracting TF-GENE groupings:
  tfList <- split(TF$FLY_TARGET_GENE,as.factor(TF$TF_SYMBOL))
  cat("| Pruning modules for genes in DroID ... ")
  L <- ontologyprune(geneList,tfList)
  # Tests:
  res <- lapply(1:length(L$geneset2),function(x){
    cat(paste0("\n|   ",names(L$geneset2)[x],": "));
    t(sapply(L$geneset1,function(y){
      buildandtest2by2(X=L$geneset2[[x]],Y=y,ALL=L$background)
    }))
  })
  names(res) <- names(tfList)
  cat("\n === DONE! ===\n\n")
  p <- sapply(res,function(x) unlist(x[,"p"]))
  OR <- sapply(res,function(x) unlist(x[,"OR"]))
  low <- sapply(res,function(x) unlist(x[,"low"]))
  high <- sapply(res,function(x) unlist(x[,"high"]))
  X11 <- sapply(res,function(x) unlist(x[,"x11"]))
  X12 <- sapply(res,function(x) unlist(x[,"x12"]))
  X21 <- sapply(res,function(x) unlist(x[,"x21"]))
  X22 <- sapply(res,function(x) unlist(x[,"x22"]))
  res <- list(p=p,OR=OR,low=low,high=high,Nmod_test=L$N1_end,Ntf_test=L$N2_end,x11=X11,x12=X12,x21=X21,x22=X22)
  return(res)
}

