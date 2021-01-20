
### DAVID Functions:
# Change-log:
# V1.0	20-06-2014	Erik	Initiation of this file
# V2.0	17-07-2014	Erik	also accepting gene lists in "DAVID"
# V2.1	18-07-2014	Erik	debug "DAVID"
# V2.2	19-07-2014	Erik	also accepting a vector of probe_names in "get.mod.members" for dat & adjust "DAVID"
# V2.3  08-05-2015	Erik	Added "sum.DAVID" to have a quick look at the results

## Defining a SubFunction for extracting the probe names of modules:
get.mod.members <- function(net,dat){
	require("WGCNA")
	require("Biobase")
	require("BiocGenerics")
	require("parallel")
	IND <- split(1:length(net$colors),as.factor(labels2colors(net$colors)))
	for(i in 1:length(IND)){
		if(is.matrix(dat)){
			names(IND[[i]]) <- rownames(dat)[IND[[i]]]
		} else {
			names(IND[[i]]) <- dat[IND[[i]]]
		}
	}
	return(IND)
}

### Defining a SubFunction for initiating my DAVID analysis:
init.DAVID <- function(){
	library("RDAVIDWebService")
	Dummi <- try(DAVIDWebService$new(email="e.b.van_den_akker@lumc.nl"),silent=TRUE) # This always gives an error the first time ...
	david <- DAVIDWebService$new(email="e.b.van_den_akker@lumc.nl")
	return(david)
}

### Defining a SubFunction for adding a gene list to DAVID WebService:
add.foreground.list <- function(con,geneIDs,listName,idType="AFFYMETRIX_3PRIME_IVT_ID"){
	# Initiate Result:
	Result <- list(succes.upload=FALSE,inDavid=NA,unmappedIds=NA,error.upload=NA,geneIDs=NA)
	# Try to add a list:
	cat(paste0("Trying to upload: \'",listName,"\' ... "))
	Dummi <- try(addList(con,geneIDs,listName=listName,idType=idType,listType="Gene"),silent=TRUE)
	# Check for error:
	if(is(Dummi,"try-error")){
		cat("FAILED! \n")
		Result$error.upload <- Dummi
	} else {
		cat("OK! \n")
		Result$succes.upload  <- TRUE
		Result$inDavid <- Dummi$inDavid
		Result$unmappedIds <- Dummi$unmappedIds
		Result$geneIDs <- geneIDs
	}
	return(Result)
}

### Defining a SubFunction for adding a list of background genes to DAVID WebService:
add.background.list <- function(con,geneIDs,listName="Background",idType="AFFYMETRIX_3PRIME_IVT_ID"){
	# Initiate Result:
	Result <- list(succes=FALSE,inDavid=NA,unMapped=NA,error=NA)
	# Try to add a list:
	cat(paste0("Trying to upload: \'",listName,"\' ... "))
	Dummi <- try(addList(con,geneIDs,listName=listName,idType=idType,listType="Background"),silent=TRUE)
	# Check for error:
	if(is(Dummi,"try-error")){
		cat("FAILED! \n")
		Result$error <- Dummi
	} else {
		cat("OK! \n")
		Result$succes  <- TRUE
		Result$inDavid <- Dummi$inDavid
		Result$unmappedIds <- Dummi$unmappedIds
	}
	return(Result)
}

### Defining a SubFunction for analyzing set listName in DAVID:
analyse.list <- function(con,listName){
	# Initiate Result:
	Result <- list(succes.analyse=FALSE,error.analyse=NA,enr=NA)
	# Check name:
	NAMES <- getGeneListNames(con)
	if(!(listName %in% NAMES)){
		Result$error <- "Genelist ot correctly uploaded"
		return(Result)
	}
	# Set list:
	setCurrentGeneListPosition(con,which(NAMES %in% listName))
	cat(paste0("Trying to analyse: \'",listName,"\' ... "))
	pathTEMP <- tempfile()
	Dummi <- try(getFunctionalAnnotationChartFile(con,fileName=pathTEMP),TRUE)
	if(is(Dummi,"try-error")){
		cat("FAILED! \n")
		Result$error.analyse <- Dummi
	} else {
		cat("OK! \n")
		Result$succes.analyse  <- TRUE
		# Read results & sort:
		Dummi <- read.delim(pathTEMP,stringsAsFactors=FALSE,header=TRUE)
		Dummi <- Dummi[sort(Dummi$Benjamini,index.return=TRUE)$ix,]
		Result$enr <- Dummi
	}
	return(Result)
}

### Defining a SubFunction for performing an DAVID enrichment analysis:
DAVID <- function(net,dat,idType="AFFYMETRIX_3PRIME_IVT_ID",annot=NULL){
	## Network results in WGCNA are not returned in its own object :(
	# Hence we need some hacking:
	if(is.list(net)){
		req.list.names <- c("colors","unmergedColors","MEs","goodSamples","goodGenes","dendrograms","TOMFiles","blockGenes","blocks","MEsOK")
		if(all(req.list.names %in% names(net))){
			# These are probably network results, thus extract module names and genes:
			cat("| Input: network \n")
			geneList  <- sapply(get.mod.members(net,dat),names)
		} else {
			# This is probably a list of genes:
			cat("| Input: list of genes \n")
			geneList <- net
		}
	}
	## Remove possible NA's:
	geneList <- sapply(geneList,function(x) x[which(!is.na(x))])
	
	## Initiate:
	david <- init.DAVID()
	david
	## Set Annotation categories:
	if(!is.null(annot)){
		NAMES <- getAllAnnotationCategoryNames(david)
		if(!all(annot %in% NAMES)){
			stop("Not all supplied annotation categories were recognized ...")
		}
	}
	setAnnotationCategories(david,annot)
	## Upload geneLists:
	if(!is.list(geneList)){
		geneList <- list(geneList=geneList)
	}
	cat("=== Uploading gene sets to DAVID ===\n")
	Result1 <- lapply(1:length(geneList),function(x){
				add.foreground.list(con=david,geneIDs=geneList[[x]],listName=names(geneList)[x],idType=idType)
			})
	names(Result1) <- names(geneList)
	cat("=== DONE! ===\n")
	## Upload background:
	cat("=== Uploading background set to DAVID ===\n")
	if(is.matrix(dat)){
		allGenes <- rownames(dat)
	} else {
		allGenes <- dat
	}
	Dummi <- add.background.list(con=david,geneIDs=allGenes,idType=idType)
	cat("=== DONE! ===\n")
	## Analyzing lists:
	cat("=== Analyzing succesfully uploaded lists ===\n")
	NAMES <- getGeneListNames(david)
	Result2 <- lapply(1:length(NAMES),function(x){
				analyse.list(con=david,listName=NAMES[x])
			})
	names(Result2) <- geneList
	Result <- lapply(1:length(geneList),function(x) c(Result1[[x]],Result2[[x]]))
	names(Result) <- names(geneList)
	cat("=== DONE! ===\n")
	return(Result)
}

### Defining a SubFunction for writing DAVID output:
write.DAVID <- function(david_res,folderOUT){
	## Create folder if not already present:
	if(!file.exists(folderOUT)){
		Dummi <- dir.create(folderOUT)
	}
	## Plot a csv per module:
	
	## Defining a helper function:
	rewrite.result <- function(i){
		mapped <- setdiff(david_res[[i]]$geneIDs,david_res[[i]]$unmappedIds)
		mapped <- data.frame(mapped,i,names(david_res)[i],stringsAsFactors=FALSE)
		result <- david_res[[i]]$enr
		## Make columns even long:
		L <- nrow(mapped)-nrow(result)
		if(L>0){
			Dummi <- data.frame(matrix(data="",ncol=ncol(result),nrow=abs(L)),stringsAsFactors=FALSE)
			colnames(Dummi) <- colnames(result)
			result <- rbind(result,Dummi)
		}
		if(L<0){
			Dummi <- data.frame(matrix(data="",ncol=ncol(mapped),nrow=abs(L)),stringsAsFactors=FALSE)
			colnames(Dummi) <- colnames(mapped)
			mapped <- rbind(mapped,Dummi)
		}
		result <- apply(cbind(mapped,result),1,function(x) paste(x,collapse="\t"))
		result <- c(paste(c("Affy","Module","Color","Category","Term","Count","%","PValue","Genes","List Total",
							"Pop Hits","Pop Total","Fold Enrichment","Bonferroni","Benjamini","FDR"),collapse="\t"),result)
		return(result)
	}
	
	pathsOUT <- paste0(folderOUT,"/",names(david_res),".txt")
	Dummi <- sapply(1:length(pathsOUT),function(x) writeLines(rewrite.result(x),pathsOUT[x]))
	return(TRUE)
}

### Defining a SubFunction for getting a summary of DAVID results:
sum.DAVID <- function(david_res){
	return(lapply(david_res,function(x) head(x$enr)[,c("Term","Fold.Enrichment","PValue","Benjamini")]))
}