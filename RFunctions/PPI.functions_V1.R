
### Setting a class to hold Protein-Protein Interactions:

# A PPIdb is a class holding interactions between genes (links), which can be scored (but not neccesarily)
setClassUnion("matrixORNULL", c("matrix", "NULL"))
setClassUnion("dataframeORNULL", c("data.frame", "NULL"))
setClass("PPIdb",representation(links="character",scores="matrixORNULL",
				metadata="data.frame"))

# An undirectedPPIdb directly inherits from PPIdb and has the additional characteristic that all
# stored links are symetric (all physical PPI's):
setClass("undirectedPPIdb",contains="PPIdb")

# A directedPPIdb directly inherits from PPIdb and has the additional characteristic that all
# stored links are directional (e.g. kinase interactions):
setClass("directedPPIdb",contains="PPIdb")

### Defining validity functions for PPIdb:
# Validity function for slot links:
.valid.ppidb.links <- function(links){
	if(!is.vector(links)){
		stop("Slot 'links' of object 'PPIdb' should be a vector")
	}
	if(!all(is.character(links))){
		stop("All entries in slot 'links' of object 'PPIdb' should be characters")
	}
	# Protein pairs making up a link are separated by a '|'. Therefore each entry should contain exactly one "|":
	if(any(!grepl("|",links,fixed=TRUE))){
		stop("All entries in slot 'links' of object 'PPIdb' should contain exactly one '|'")
	}
	return(TRUE)
}
# Validity function for slot scores:
.valid.ppidb.scores <- function(scores){
	if(!is(scores,"matrixORNULL")){
		stop("Slot 'scores' of object 'PPIdb' should either be a matrix or NULL")
	}
	if(is.matrix(scores)){
		if(!all(is.numeric(scores))){
			stop("All entries in slot 'scores' of object 'PPIdb' should be numeric")
		}
	}
	return(TRUE)
}
# Validity function for slot metadata:
.valid.ppidb.metadata <- function(metadata){
	if(!is.data.frame(metadata)){
		stop("Slot 'metadata' of object 'PPIdb' should be a supplied as a 'data.frame'")
	}
	if(ncol(metadata)!=2){
		stop("Slot 'metadata' of object 'PPIdb' should be a 'data.frame' of two columns")
	}
	if(!all(colnames(metadata)==c("name","value"))){
		stop("Slot 'metadata' of object 'PPIdb' should be a 'data.frame' with header: 'name','value'")
	}
}
# Validity function for PPIdb object:
.valid.PPIdb <- function(object){
	# 1) Check slots of PPIdb object:
	# 1a) links:
	Dummi <- .valid.ppidb.links(object@links)
	# 1b) scores:
	Dummi <- .valid.ppidb.scores(object@scores)
	# 1c) metadata:
	Dummi <- .valid.ppidb.metadata(object@metadata)
	# 2) Concistency between slots:
	# 2a) matching dimensions:
	if(!is.null(object@scores)){
		if(length(object@links)!=nrow(object@scores)){
			stop("Number of entries in 'scores' should match 'links'")
		}
	}
}
# Set validity function for PPIdb:
setValidity("PPIdb",.valid.PPIdb)

### Low level accessors:

### Defining a standardGeneric for accessing links of an object:
setGeneric("links",function(object,...) standardGeneric("links"))

## Defining a Method for accessing links of a PPIdb object:
setMethod("links",signature(object="PPIdb"),
		function(object,...){
			Dummi <- .valid.ppidb.links(object@links)
			if(!is.null(object@links)){
				Dummi <- t(simplify2array(strsplit(object@links,"|",fixed=TRUE)))
				colnames(Dummi) <- c("protID1","protID2")
				return(Dummi)
			} else {
				return(NULL)
			}
		}
)

### Defining a standardGeneric for accessing scores of an object:
setGeneric("scores",function(object,...) standardGeneric("scores"))

## Defining a Method for accessing scores of a PPIdb object:
setMethod("scores",signature(object="PPIdb"),
		function(object,...){
			Dummi <- .valid.ppidb.scores(object@scores)
			Dummi <- object@scores
			if(is.null(Dummi)){
				Dummi <- matrix(nrow=length(object@links),ncol=0)
			}
			rownames(Dummi) <- object@links
			return(Dummi)
		}
)

### Defining a standardGeneric for accessing metadata of an object:
setGeneric("meta",function(object,...) standardGeneric("meta"))

## Defining a Method for accessing metadata of a PPIdb object:
setMethod("meta",signature(object="PPIdb"),
		function(object,...){
			Dummi <- .valid.ppidb.metadata(object@metadata)
			Dummi <- object@metadata
			return(Dummi)
		}
)

### Defining a low level constructor for PPIdb:
.PPIdb <- function(links,scores,metadata,directed){	
	if(directed){
		Result <- new("directedPPIdb",links=links ,scores=scores, metadata = metadata)
	} else {
		Result <- new("undirectedPPIdb",links=links ,scores=scores, metadata = metadata)
	}
	return(Result)
}

### Defining a function for making a PPIdb:
makePPIdb <- function(ppi.df,meta.df=NULL,directed=FALSE){
	### Checking input:
	# 1) Data.frame with ppi information:
	# 1a) PPI data:
	if(!is.data.frame(ppi.df)){
		stop("Please supply a 'data.frame' for 'ppi.df'")
	}
	.req.colnames.ppi.df <- c("protID1","protID2")
	if(all(!(colnames(ppi.df) %in% .req.colnames.ppi.df))){
		stop("Please supply a 'data.frame' for 'ppi.df' with header containing 'protID1' and 'protID2'")
	}
	if(all(!is.character(ppi.df$protID1))){
		stop("Please supply character values for slot 'protID1' of 'ppi.df'")
	}
	if(all(!is.character(ppi.df$protID2))){
		stop("Please supply character values for slot 'protID2' of 'ppi.df'")
	}
	if(ncol(ppi.df)>2){
		# Scores are provided:
		SCORES <- ppi.df[,setdiff(colnames(ppi.df),.req.colnames.ppi.df)]
		for(i in 1:ncol(SCORES)){
			if(all(is.logical(SCORES[,i]))){
				SCORES[,i] <- as.numeric(SCORES[,i])
			}
			if(!all(is.numeric(SCORES[,i]))){
				stop("Please supply numerical values for scores of ppi's in 'ppi.df'")
			}
		}
		SCORES <- as.matrix(SCORES)
	} else {
		SCORES <- NULL
	}
	# 1b) metadata:
	if(!is(meta.df,"dataframeORNULL")){
		stop("Please supply a either NULL or a 'data.frame' for 'meta.df'")
	}
	if(is.data.frame(meta.df)){
		.req.colnames.meta.df <- c("name","value")
		if(all(!(colnames(meta.df) %in% .req.colnames.meta.df))){
			stop("Please supply a 'data.frame' for 'meta.df' with header containing 'name' and 'value'")
		}
	}
	# 1c) directed:
	if(!is.logical(directed)){
		directed <- FALSE
	}
	### Rewrite links:
	if(!directed){
		L <- apply(ppi.df[,c("protID1","protID2")],1,function(x) paste(sort(x),collapse="|"))
	} else {
		L <- paste(ppi.df$protID1,ppi.df$protID2,sep="|")
	}
	### Creating minimal metadata:
	auto_metadata <- matrix(c(
					"N_Links"			, length(L),
					"N_Nodes"  			, length(unique(c(ppi.df$protID1,ppi.df$protID2))),
					"N_Scores"			, ifelse(is.null(SCORES),0,ncol(SCORES)),
					"dbCreatedBy"		, "PPI.functions",
					"Creation time" 	, gsub(" ","_",as.character(Sys.time()))),ncol=2,byrow=TRUE)
	colnames(auto_metadata) <- c("name", "value")
	auto_metadata <- as.data.frame(auto_metadata,stringsAsFactors=FALSE)
	if(is.null(meta.df)){
		metadata <- auto_metadata
	} else {
		meta.df <- meta.df[,c("name","value")]
		auto_metadata <- auto_metadata[which(!(auto_metadata$name %in% meta.df$name)),]
		metadata <- rbind(meta.df,auto_metadata)
	}
	### Creating:
	Result <- .PPIdb(links = L, scores = SCORES, metadata = metadata, directed = directed)
	return(Result)
}

## Setting method show:
setMethod("show", "PPIdb",
		function(object){
			cat("--- Protein-Protein Interaction database: ---\n")
			# 1) show metadata:
			if(is(object,"directedPPIdb")){
				cat("| directed: TRUE \n")
			}
			if(is(object,"undirectedPPIdb")){
				cat("| directed: FALSE \n")
			}
			for(i in seq_len(nrow(object@metadata))){
				cat("| ", object@metadata[i, "name"], ": ", 
						object@metadata[i, "value"],"\n", sep="")
			}
			cat("| \n")
			# 2) show example output:
			# 2a) header:
			cat(paste("| First few lines: \n"))
			cat("| \n")
			# 2b) entries:
			maxdim <- ifelse(length(object@links)>3,3,length(object@links))
			if(maxdim!=0){
				Example <- scores(object[1:maxdim,,drop=FALSE])
				rownames(Example) <- paste0("| ",rownames(Example))
				show(Example)
			}
		}
)

## Setting method [:
setMethod("[", "PPIdb",
		function(x,i,drop=F) {
			if(is.character(i)){
				i <- match(i,x@links)
			}
			x@links  <- x@links[i,drop=F]
			x@scores <- x@scores[i,,drop=F]
			if(all(i==0)){
				x@metadata$value[which(x@metadata$name=="N_Links")] 	<- 
						as.character(0)
			} else {
				x@metadata$value[which(x@metadata$name=="N_Links")] 	<- 
						as.character(length(i))
			}
			x@metadata$value[which(x@metadata$name=="N_Nodes")] 	<-
						length(unique(unlist(strsplit(x@links,split="|"))))
			return(x)
		}
)

## Setting method dim:
setMethod("dim",signature(x="PPIdb"),
		function(x){
			return(dim(x@scores))
		}
)

### Defining a Subfunction for mapping affyIDs to fbIDs (FlyBaseIDs). 
# Only uniquely mapped probes are returned!
# Type rev=TRUE, for an inverse mapping.
Affy2FlyBase <- function(affyIDs=NULL,rev=FALSE){
	require("drosophila2.db")
	if(is.null(affyIDs)){
		affyIDs <- grep("AFFX",keys(drosophila2.db),value=TRUE,invert=TRUE)
	}
	suppressWarnings(annot <- select(drosophila2.db,keys=affyIDs, 
					cols=c("FLYBASE"), keytype="PROBEID"))
	# Search for duplicate mappings & prune:
	Dummi <- rle(sort(annot[,1]))
	dupID <- Dummi$values[which(Dummi$lengths>1)]
	if(length(dupID)>0){
		cat("Removing ",length(dupID)," probes mapping to multiple genes \n")
	}
	annot <- annot[which(!(annot$PROBEID %in% dupID)),]
	# Make mapping:
	if(nrow(annot)>1){
		if(!rev){
			# Make a forward mapping: values are fbIDs & names are affyIDs:
			annot2 <- annot$FLYBASE
			names(annot2) <- annot$PROBEID
		} else {
			# Make a reverse mapping: values are affyIDs & names are fbIDs:
			annot2 <- annot$PROBEID
			names(annot2) <- annot$FLYBASE
		}
		return(annot2)
	} else {
		return(NULL)
	}	
}

### Defining a SubFunction for selecting links for a set of genes:
getModulePPIs <- function(geneList,ppidb,simplify=TRUE){
	if(!is.list(geneList)){
		geneList <- list(list1=geneList)
	}
	L <- links(ppidb)
	IND <- lapply(geneList,function(x) which((x %in% L[,1]) & (x %in% L[,2])))
	for(i in 1:length(IND)){
		names(IND[[i]]) <- ppidb@links[IND[[i]]]
	}
	if(length(IND)==1&simplify){
		return(IND[[1]])
	}
	return(IND)
}

### Defining a SubFunction for perfomring a PPI enrichment test:
ppi.enr <- function(geneList,ppidb,allGenes=NULL,Ntimes=10000,seed=12345){
	require("plyr")
	cat(paste0("=== Performing PPI Enrichment analyses on ",length(geneList)," gene sets ===\n"))
	if(!is.list(geneList)){
		geneList <- list(list1=geneList)
	}
	if(is.null(allGenes)){
		allGenes <- unlist(geneList)
	}
	## Compute observed:
	Ngenes <- sapply(geneList,length)
	Obs_IND <- getModulePPIs(geneList,ppidb,simplify=FALSE)
	if(is.null(scores(ppidb))){
		SCORES <- cbind(scores(ppidb),rep(1,length(ppidb@links)))
		colnames(SCORES)[1] <- "links"
	} else {
		SCORES <- scores(ppidb)
	}
	Obs <- t(sapply(Obs_IND,function(x) colSums(SCORES[x,])))
	
	## Compute background distributions for each unique size in geneList:
	if(!is.null(seed)){
		set.seed(seed)
	}
	## Sample genes:
	cat(paste0("| Computing background distributions: \n| "))
	gene_IND <- replicate(Ntimes,sample(1:length(allGenes),max(Ngenes)))
	## Find link indices:
	Ugenes <- unique(Ngenes)
	Exp_IND <- lapply(Ugenes,function(x){
		cat("*");getModulePPIs(alply(gene_IND[c(1:x),],2,function(y) allGenes[y]),ppidb,simplify=FALSE)})	
	names(Exp_IND) <- Ugenes
	cat(" Done!\n")
	
	## Defining a SubFunction for looking up scores in equally sized lists of link indices:
	lookup.list.equal <- function(equalsized_lists,scores){
		Nlinks <- length(equalsized_lists[[1]])
		if(Nlinks==0){
			SumS <- matrix(0,ncol=ncol(scores),nrow=length(equalsized_lists))
			colnames(SumS) <- colnames(scores)
		} else {
			## Look up as a vector:
			S <- scores[unlist(equalsized_lists,use.names=FALSE),,drop=FALSE]
			## Rebuild matrices and compute sum.statistic:
			SumS <- apply(S,2,function(x) colSums(matrix(x,nrow=Nlinks,byrow=FALSE)))
		}
		return(SumS)
	}

	## Defining a SubFunction for efficient lookup in ppidb:
	lookup.list <- function(linkIND,scores){
		## Divide links in groups of links who are equally sized:
		L <- sapply(linkIND,length)
		IND_IND <- split(L,as.factor(L))
		## Go over these batches of equally sized link groups and compute their summed statistics per batch: 
		SumS <- do.call("rbind",lapply(IND_IND,function(x) lookup.list.equal(linkIND[names(x)],scores=scores)))
		rownames(SumS) <- NULL
		return(SumS)
	}
	
	## Defining a SubFunction for getting the statistics:
	gimme.pees <- function(foreground,background){
		## Compute Z scores:
		M  <- colMeans(background)
		SD <- apply(background,2,sd)
		Z  <- as.vector((foreground-M)/SD)
		p.Z <- pnorm(-abs(Z))
		## Compute permutation.p:
		Nbeat_l <- sapply(1:ncol(foreground),function(x) length(which(background[,x]<=foreground[,x])))
		Nbeat_h <- sapply(1:ncol(foreground),function(x) length(which(background[,x]>=foreground[,x])))
		pval_l  <- (Nbeat_l+1)/(Ntimes+1)
		pval_h  <- (Nbeat_h+1)/(Ntimes+1)
		Result <- data.frame(obs=as.vector(foreground),Z.obs=Z,M.back=M,SD.back=SD,p.obs=p.Z,
					Nbeat_low=Nbeat_l,pval_low=pval_l,Nbeat_high=Nbeat_h,pval_high=pval_h)
		return(Result)
	}
	
	## Lookup scores and compute background distribution of sum statistics:
	SumS_Exp <- lapply(Exp_IND,function(x) lookup.list(x,scores=SCORES))
	## Go over geneList and compute tests:
	Result <- lapply(1:length(geneList),function(x){ 
			gimme.pees(Obs[x,,drop=FALSE],SumS_Exp[[as.character(length(geneList[[x]]))]])})
	names(Result) <- names(geneList)
	return(Result)
}



