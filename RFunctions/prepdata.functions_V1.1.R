
### Prepare Data Functions:

### Aim: To collect all functions relating to preprocessing expression data 
# (limited to fly Affy arrays for the moment).

# Change-log:
# V1.0	07-05-2015	Erik	Initiation of this file

### Defining a SubFunction for matching pheno to the data:
match.pheno <- function(pheno,pathsIN_CEL){
    cat("=== Matching phenotype matrix to CEL files === \n")
    pheno <- pheno[sort(pheno$sampleID,index.return=T)$ix,]
    pheno$Order.pheno <- 1:nrow(pheno)
    CEL_ID <- data.frame(   pathIN=pathsIN_CEL,
                            sampleID=gsub(".CEL","",sapply(strsplit(basename(pathsIN_CEL),"_"),function(x) x[1])),
                            Order.cell=1:length(pathsIN_CEL),    
                    		stringsAsFactors=FALSE)
    res <- merge(pheno,CEL_ID,by.x="sampleID",by.y="sampleID")
    rownames(res) <- res$sampleID
    cat(paste0("| # Entries in pheno matched: ",nrow(res)," out of ",nrow(pheno),"\n"))
    cat(paste0("| # Supplied CELL files matched: ",nrow(res)," out of ",length(pathsIN_CEL),"\n"))
    return(res)
}

### Defining a SubFunction for preprocessing Flyatlas data:
prep.fa.pheno <- function(pathIN_Pheno){
    cat("=== Preparing Flyatlas phenotypes === \n")
    ### 1) Reading data:
    ## Phenotypes:
    pheno <- readLines(pathIN_Pheno)
    # Remove original samples of those who are redone:
    pheno <- pheno[grep("Adult Female Spermatheca Mated rep",pheno,invert=TRUE)]
    # Get ID:
    sampID <- sapply(strsplit(pheno,"\t"),function(x) x[1])
    # Origin:
    origin  <- rep("",length(pheno))
    origin[grep("Adult",pheno)]             <- "ADULT"
    origin[grep("Larval|Larvae",pheno)]     <- "LARVAL"
    origin[grep("Cells",pheno)]             <- "CELL"
    origin <-factor(origin)
    # Replicate:
    Dummi <- sapply(strsplit(pheno,"\t"),function(x) x[2])
    rep <- substr(Dummi,nchar(Dummi)-3,nchar(Dummi))
    # Tissue:
    tissue <- gsub(" ","_",gsub(" [(]REDONE[)] biological rep[1-4]| rep[1-4]| biological rep[1-4]","",Dummi))
    tissue <-factor(tissue)
    pheno <- data.frame(sampleID=sampID,origin=origin,tissue=tissue,rep=rep,stringsAsFactors=FALSE)
    cat(paste0("| # Entries in pheno read: ",nrow(pheno),"\n"))
    return(pheno[,c("sampleID","origin","tissue","rep")])
}

### Defining a SubFunction for to read, normalize and filter data:
norm.dat <- function(pheno,pathsIN_CEL){
	require("Biobase")
	require("gcrma")
	require("affy")
	require("genefilter")
	## 1] Matching Pheno to CEL files:
	pheno <- match.pheno(pheno=pheno,pathsIN_CEL=pathsIN_CEL)
	## 2] Read Expression Data:
	COLNAMES <- setdiff(colnames(pheno),c("Order.pheno","Order.cell","pathIN"))
	abatch <- read.affybatch(filenames=pheno$pathIN,phenoData=as(pheno[,COLNAMES],"AnnotatedDataFrame"))
	# 3] Normalize:
	eset <- gcrma(abatch)
	# 4] Remove "^AFFX" probes:
	Expr  <- nsFilter(eset, require.entrez = FALSE, remove.dupEntrez = FALSE, var.func = IQR, 
		var.cutoff = 0.5, var.filter = FALSE, filterByQuantile=TRUE,feature.exclude="^AFFX")$eset
	return(Expr=Expr)
}

