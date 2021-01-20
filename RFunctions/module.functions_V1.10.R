
### Module Functions:
# Change-log:
# V1.0	18-07-2014	Erik	Initiation of this file
# V1.1	19-07-2014	Erik	Pimping "plot.module.size", add "plot.module.overlap"
# V1.2  23-07-2014	Erik	debug introduced incompatibility (V1.0->V1.1) in "plot.module.overlap"
# V1.3	19-01-2015	Erik	Add functions: "get.module.me", "get.sample.similarity","plot.sample.similarity.vs.module.eg","plot.module.heatmap"
# V1.4	21-01-2015	Erik	debug introduced incompatibility (V1.2->V1.3) in "plot.sample.dendro"
# V1.5	03-02-2015	Erik	Add: 'rec.eliminate.sample'
# V1.6  20-02-2015	Erik	Add: samplenames to ME matrix obtained with 'get.module.eg'
#							            Add: Function for computing the median of a module: 'get.module.mg'
# V1.6	07-05-2015	Erik	Add: Function for performing tissue enrichments on expression abundance: 'mod.tis.enr'
#							            Add: Function for plotting tissue enrichments: 'plot.tis.enr'
# V1.7	21-05-2015	Erik	Fixed bug in: "plot.sample.dendro": Outliers were not detected correctly!
#							            Add: Function for listing outliers per module: 'get.module.outliers'
# V1.8  17-07-2015  Erik  Fixed bug in: "mod.tis.enr" added pseudo-counts to prevent computed odd's ratio's to go to infinity!
#							            This makes the test more conservative, but also capable of handling perfect overlaps.
# V1.9  31-07-2015  Erik  Solved small problem with flashClust dependencies in R 3.2.1
# V1.10 01-08-2015  Erik  Add options to change color and shape of points used in "plot.sample.similarity.vs.module.eg"
#                         Add options to make a selection of modules plotted

## Defining a SubFunction for obtaining a sample similarity measure:
get.sample.similarity <- function(dat){
	require("WGCNA")
	## Compute Sample Similarities:
	A <- adjacency(dat,type="distance")
	Z.k <- scale(as.numeric(apply(A,2,sum))-1)
	return(Z.k)
}

## Defining a Subfunction for plotting sample similarities:
plot.sample.similarity <- function(dat,Zcut=2.5,main=NULL){
	## Compute Sample Similarities:
	Z.k <- get.sample.similarity(dat)
	## Settings Figure:
	XLIM <- c(0,1)
	YRANGE <- pretty(Z.k)
	YLIM   <- c(-max(abs(YRANGE),(Zcut+0.5)),max(abs(YRANGE),(Zcut+0.5)))
	YTICK  <- pretty(YLIM)
	plot.new()
	plot.window(xlim=XLIM,ylim=YLIM)
	## Indicate Outliers:
	Ind0 <- which((abs(Z.k) > Zcut))
	Ind1 <- which(!(abs(Z.k) > Zcut))
	## Draw Boundaries:
	lines(c(0,1),c(0,0),lwd=2,col="black")
	lines(c(0,1),c(-Zcut,-Zcut),lty=2,lwd=2,col="black")
	lines(c(0,1),c(Zcut,Zcut),lty=2,lwd=2,col="black")
	## Draw samples:
	if(length(Ind0)>0){
		points(x=rep(0.5,length(Ind0)),y=Z.k[Ind0],bg=rgb(1,0,0,0.5),pch=21)
		text(x=0.6,y=Z.k[Ind0],colnames(dat)[Ind0],col=rgb(1,0,0),font=2,pos=4)
		text(x=0.4,y=Z.k[Ind0],paste0("Zk=",signif(Z.k,3)[Ind0]),col=rgb(1,0,0),font=2,pos=2)
	}
	if(length(Ind1)>0){
		points(x=rep(0.5,length(Ind1)),y=Z.k[Ind1],bg=rgb(0,0,1,0.5),pch=21)
	}
	## Do some Figure Pimping:
	legend("topleft",legend=c("Outlier: |Zk|>Zcut","Normal: -Zcut<=Zk<=Zcut"),
			fill=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)))
	axis(2,at=YTICK,labels=as.character(YTICK))
	box()
	if(is.null(main)){
		main <- "Sample Similarity based on Euclidian Distance"
	}
	title(main=main,ylab="Z-Scaled Similarity")
	## Return names of outliers if any present:
	if(length(Ind0)>0){
		return(colnames(dat)[Ind0])
	} else {
		return(NULL)
	}
}

## Defining a SubFunction for plotting sample dendrograms:
plot.sample.dendro <- function(dat,Zcut=2.5,pheno=NULL,main=NULL){
	require("WGCNA")
	require("flashClust")
	## Compute Sample Similarities on basis of expression data & determine outliers:
	A <- adjacency(dat,type="distance")
	Z.k <- scale(as.numeric(apply(A,2,sum))-1)
	outlierColor=ifelse(abs(Z.k)>Zcut,"red","black")
	## Compute clustering:
	tree = flashClust(as.dist(1-A), method = "average")
	if(!is.null(pheno)){
		traitColors=data.frame(labels2colors(pheno))
		dimnames(traitColors)[[2]]=names(pheno)
		datColors=data.frame(outlier=outlierColor,traitColors)
		if(is.null(main)){
			main <- "Sample dendrogram and trait heatmap"
		}
	} else {
		datColors=data.frame(outlier=outlierColor)
		if(is.null(main)){
			main <- "Sample dendrogram"
		}
	}
	# Plot the sample dendrogram and the colors underneath.
	plotDendroAndColors(tree,groupLabels=names(datColors),
	colors=datColors,main=main)
}

## Defining a Subfunction for computing modules:
comp.modules <- function(dat,type="Tina"){
	if(type=="Tina"){
		net <- blockwiseModules(t(dat),
				maxBlockSize=nrow(dat)+1,
				corType="pearson",
				networkType="signed",
				power=12,
				minModuleSize=30,
				mergeCutHeight=0.25,
				numericLabels=TRUE,
				saveTOMs=FALSE,
				pamRespectsDendro=FALSE)
	}
	if(type=="Erik"){
		net <- blockwiseModules(t(dat), 
				maxBlockSize=nrow(dat)+1,
				corType="bicor",
				networkType="signed",
				power=12,
				numericLabels=TRUE)
	}
	return(net)
}
		
## Defining a SubFunction for plotting modules:
plot.module.dendro <- function(net){
	modCOL <- labels2colors(net$colors)
	plotDendroAndColors(net$dendrograms[[1]],
		colors=data.frame(modCOL),
		groupLabels=c("Module \n Membership"),
		dendroLabels=FALSE,
		hang=0.03,
		addGuide=TRUE,
		guideHang=0.05)
}

## Defining a SubFunction for extracting the module eigengenes:
get.module.eg <- function(dat,net){
	ME <- orderMEs(moduleEigengenes(t(dat),colors=labels2colors(net$colors))$eigengenes)
	rownames(ME) <- colnames(dat)
	return(ME)
}

## Defining a SubFunction for extracting the module mediangenes:
get.module.mg <- function(dat,net){
	IND <- split(1:length(net$colors),as.factor(labels2colors(net$colors)))
	MM <- sapply(IND,function(x) apply(dat[x,],2,median))
	return(MM)
}

## Defining a SubFunction to plot a dendrogram of module eigengenes:
plot.ME.dendro <- function(dat,net){
	## Compute and order module eigengenes:
	#ME <- orderMEs(net$MEs)
	ME <- get.module.eg(dat,net)
	## Plot:
	plotEigengeneNetworks(ME,"",marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)
}

## Defining a SubFunction for a barplot to vizualize sizes of genes per module:
plot.module.size <- function(net,main=NULL){
	Dummi <- sort(sapply(split(net$colors,as.factor(labels2colors(net$colors))),length),decreasing=TRUE)
	XRANGE <- pretty(Dummi)
	XLIM <- c(min(XRANGE),max(XRANGE))
	if(is.null(main)){
		main <- "Sizes of Modules"
	}
	par(oma=c(5.0,5.0,4.0,2.0))
	bp <- barplot(Dummi,col=names(Dummi),horiz=TRUE,las=1,main=main,xlab="# of probes",xlim=XLIM)
	text(Dummi, unlist(bp), Dummi, pos=4, col="red",cex=0.8,font=2)
}

## Defining a SubFunction for vizualizing module overlap: 
plot.module.overlap <- function(net1,net2,main=NULL){
	# Compute numbers of intersecting probes:
	IND_1 <- split(1:length(net1$colors),as.factor(labels2colors(net1$colors)))
	IND_2 <- split(1:length(net2$colors),as.factor(labels2colors(net2$colors)))
	modOverlap <- sapply(IND_1,function(x) sapply(IND_2, function(y) length(intersect(x,y))))
	tot1 <- as.vector(sort(unlist(IND_1)))
	tot2 <- as.vector(sort(unlist(IND_2)))
	
	## Defining a SubFunction for computing a fisher's test:
	fish <- function(x1,x2,tot1,tot2){
		tot <- intersect(tot1,tot2)
		x1 <- x1[which(x1 %in% tot)]
		x2 <- x2[which(x2 %in% tot)]
		A <- intersect(x1,x2)
		B <- intersect(x1,setdiff(tot,x2))
		C <- intersect(x2,setdiff(tot,x1))
		D <- setdiff(tot,intersect(x1,x2))
		return(fisher.test(matrix(c(length(A),length(B),length(C),length(D)),nrow=2,byrow=T))$p.value)
	}
	# Compute significance of overlap:
	modPVal    <- sapply(IND_1,function(x) sapply(IND_2, function(y) fish(x1=x,x2=y,tot1,tot2)))
	# Preparing output:
	textMatrix = paste(modOverlap, "\n(",signif(modPVal, 1), ")", sep = "")
	dim(textMatrix) = dim(modPVal)
	par(mar = c(6, 8.8, 3, 2.2))
	Dummi <- -log10(modPVal)
	Dummi[which(modPVal==0)] <- NA
	# Plot:
	if(is.null(main)){
		main <- "Module Overlap"
	}
	labeledHeatmap(Matrix = Dummi,
		xLabels = names(IND_1),
		yLabels = names(IND_2),
		ySymbols = names(IND_2),
		colorLabels = FALSE,
		colors=greenWhiteRed(100)[51:100], 
		textMatrix=textMatrix,setStdMargins=FALSE,cex.text=0.5,main=main)
}

### Defining a SubFunction for plotting module eigen-genes vs. Z.k:
plot.sample.similarity.vs.module.eg <- function(dat,net,Zcut=2.5,hold=TRUE,col=NULL,pch=NULL,module=NULL){
	## Defining a Subfunction for each plot:
	plot.me_zk <- function(me,zk,name,Zcut,col,pch){
		## Settings Figure:
		XRANGE <- pretty(pretty(me)*1.1)
		XLIM   <- c(-abs(XRANGE),abs(XRANGE))
		XTICK  <- pretty(XLIM)
		XLIM   <- c(min(XTICK),max(XTICK))
		YRANGE <- pretty(zk)
		YLIM   <- c(-max(abs(YRANGE),(Zcut+0.5)),max(abs(YRANGE),(Zcut+0.5)))
		YTICK  <- pretty(YLIM)
		plot.new()
		plot.window(xlim=XLIM,ylim=YLIM)
		## Indicate Outliers on X and Y:
		Ind0 <- which((as.vector(abs(scale(me)) > Zcut))|(as.vector(abs(Z.k) > Zcut)))
		Ind1 <- which(!(abs(scale(me)) > Zcut)&!(abs(Z.k) > Zcut))
		## Draw Boundaries:
		lines(XLIM,c(0,0),lwd=2,col="black")
		lines(XLIM,c(-Zcut,-Zcut),lty=2,lwd=2,col="black")
		lines(XLIM,c(Zcut,Zcut),lty=2,lwd=2,col="black")
		## Draw samples:
		if(is.null(col)){
		  col <- rgb(0,0,1,0.5)
		}
		if(is.null(pch)){
		  pch <- 21
		}
		if(length(Ind0)>0){
			DEV1 <- 0.01*(XLIM[2]-XLIM[1])
			DEV2 <- 0.05*(YLIM[2]-YLIM[1])
			points(x=me[Ind0],y=Z.k[Ind0],bg=col,pch=pch,cex=2)
			text(x=me[Ind0]+DEV1,y=Z.k[Ind0],colnames(dat)[Ind0],col=rgb(1,0,0),font=2,pos=4)
			text(x=me[Ind0]-DEV1,y=Z.k[Ind0],paste0("Zk=",signif(Z.k,3)[Ind0]),col=rgb(1,0,0),font=2,pos=2)
			text(x=me[Ind0]-DEV1,y=Z.k[Ind0]-DEV2,paste0("Zme=",signif(scale(me),3)[Ind0]),col=rgb(1,0,0),font=2,pos=2)
		}
		if(length(Ind1)>0){
			points(x=me[Ind1],y=Z.k[Ind1],bg=col,pch=pch,cex=1)
		}
		## Do some Figure Pimping:
		axis(1,at=XTICK,labels=as.character(XTICK))
		axis(2,at=YTICK,labels=as.character(YTICK))
		box()
		TEST <- cor.test(me,zk)
		title(main=paste0(name," COR=",as.character(signif(cor(me,zk),3))," p=",signif(TEST$p.value,3)),ylab="Z.k",xlab="Module Expression")	
		if(hold){
			cat ("Press [enter] to continue")
			line <- readline()
		}
	}
	# Get similarities between samples:
	Z.k <- get.sample.similarity(dat)
	# Get module eigengenes:
	ME <- get.module.eg(dat,net)
	## Make a scatterplot for each module:
	if(is.null(module)){
  	for(i in 1:ncol(ME)){
  		plot.me_zk(ME[,i],Z.k,colnames(ME)[i],Zcut,col,pch)
  	}
	} else {
	  for(i in 1:length(module)){
	    plot.me_zk(ME[,paste0("ME",module[i])],Z.k,paste0("ME",module[i]),Zcut,col,pch)  
	  }
	}
}

### Defining a SubFunction for listing module eigengene outliers:
get.module.outliers <- function(dat,net,phen,Zcut){
	# Get module eigengenes:
	ME <- get.module.eg(dat,net)
	# List outliers:
	OUT <- apply(ME,2,function(x){
		data.frame(phen,Z=scale(x),stringsAsFactors=FALSE)[which((abs(scale(x))) > Zcut),]
	})
	names(OUT) <- colnames(ME)
	return(OUT)
}

### Defining a SubFunction for plotting module heatmaps:
plot.module.heatmap <- function(dat,net,col){
	require("gplots")
	if(!(col %in% unique(labels2colors(net$colors)))){
		stop(paste0("'col' must be one of: ",paste(unique(labels2colors(net$colors)),collapse=", ")))
	}
	IND <- which(labels2colors(net$colors)==col)
	par(oma=c(0,0,0,0))
	heatmap.2(dat[IND,],Rowv=FALSE,Colv=FALSE,dendrogram="none",
			scale="row",key=TRUE,trace="none",labRow=NULL,labCol=NULL)
}

## Defining a SubFunction to recursively eliminate outlier samples:
rec.eliminate.sample <- function(dat,Zcut=2.5,main=NULL,plot=TRUE){
	cat("=== Recursively eliminating outlier samples ===\n")
	continue <- TRUE
	i <- 0
	sampREMOVE <- NULL
	while(continue==TRUE){
		dat <- dat[,which(!(colnames(dat) %in% sampREMOVE))]
		Zk <- get.sample.similarity(dat)
		sampREMOVE <- colnames(dat)[which.max(abs(Zk))]
		if(max(abs(Zk))>=Zcut){
			cat(paste0("| Eliminating [",i,"]: ",sampREMOVE," (Zk=",signif(Zk[which.max(abs(Zk))],3),")\n"))
		} else {
			cat(paste0("| Stopping elimination process at [",i,"] with: ",sampREMOVE," (Zk=",signif(Zk[which.max(abs(Zk))],3),")\n"))
			sampREMOVE <- NULL
		}
		if(plot){
			if(is.null(main)){
				new_main <- paste0("Elination Round: ",i)
			} else {
				new_main <- paste0(main," [",i,"]")
			}
			Dummi <- plot.sample.similarity(dat=dat,Zcut=Zcut,main=new_main)
		}
		if((i>0)&(length(sampREMOVE)==0)){
			continue <- FALSE
		}
		i <- i+1
	}
}

## Defining a SubFunction for computing tissue enrichments:
mod.tis.enr <- function(net,dat,pathIN_TA,Ecut=5,background="all"){

	## Defining a SubFunction for performing a 2x2 test:
	buildandtest2by2 <- function(X,Y,ALL){
		cat("*")
		if(is.null(X)){
			X <- ""
		}
		x11 <- length(intersect(X,Y))
		x12 <- length(setdiff(Y,X))+1
		x21 <- length(setdiff(X,Y))+1
		x22 <- length(setdiff(ALL,c(X,Y)))

		res <- fisher.test(matrix(c(x11,x12,x21,x22),nrow=2,byrow=TRUE))
		return(list(p=res$p.value,OR=res$estimate[[1]],low=res$conf.int[1],high=res$conf.int[2]))
	}
	require("Biobase")
	cat("=== Performing Tissue Enrichments === \n")
	## Load prepared Tissue Atlas data:
	cat(paste("| Preparing Tissue Atlas Data ... "))
	FA <- get(load(pathIN_TA))
	## Summarize by taking the median over replicates:
	pheno_FA <- as(phenoData(FA),"data.frame")
	IND <- split(1:nrow(pheno_FA),droplevels(pheno_FA$tissue))
	dat_FA <- sapply(IND,function(x) apply(exprs(FA)[,x],1,median))
	cat("OK! \n")
	## Apply cutoff:
	cut <- dat_FA>=Ecut
	N1_tis <- colSums(cut)
	## Test:
	allow.bg <- c("expressed","all")
	if(!(background %in% allow.bg)){
		stop(paste("Please select for 'background' one of: '",paste(allow.bg,collapse="','"),"'",sep=""))
	}
	if(background=="all"){
		bg <- rownames(cut)
	} else {
		keepIND <- which(rowSums(cut)>0)
		bg <- names(keepIND)
		cut <- cut[keepIND,]
	}
	cat(paste0("| Setting background : '",background,"' [N = ",length(bg),"] \n"))
	cat(paste0("| Testing for enrichments per tissue at Absolute Expression >= ",Ecut,"..."))
	geneList  <- sapply(get.mod.members(net,dat),names)
	N0_mod <- sapply(geneList,length)
	if(background=="expressed"){
		cat(paste("| Selecting genes expressed in ANY tissue of the atlas ..."))
		geneList <- sapply(geneList,function(x) x[which(x %in% bg)])
	}
	res <- lapply(1:ncol(cut),function(x){
		cat(paste0("\n|   ",colnames(cut)[x]," [# Expr >= ",Ecut,": ",length(which(cut[,x])),"]: \n"));
		t(sapply(geneList,function(y){
			buildandtest2by2(names(which(cut[,x])),y,bg)
		}))
	})
	cat("\n === DONE! ===\n\n")
	p <- sapply(res,function(x) unlist(x[,"p"]))
	colnames(p) <- colnames(cut)
	OR <- sapply(res,function(x) unlist(x[,"OR"]))
	colnames(OR) <- colnames(cut)
	low <- sapply(res,function(x) unlist(x[,"low"]))
	colnames(low) <- colnames(cut)
	high <- sapply(res,function(x) unlist(x[,"high"]))
	colnames(high) <- colnames(cut)
	N1_mod <- sapply(geneList,length)
	res <- list(p=p,OR=OR,low=low,high=high,Ngenes_input_set=N0_mod,Ngenes_test_set=N1_mod,
					Ngenes_test_tissue=N1_tis)
	return(res)
}

## Defining a SubFunction for plotting tissue enrichments:
plot.tis.enr <- function(enr,pmin=1e-15){
	require("gplots")
	pmat <- enr$p
	pmat[which(pmat<=pmin)] <- pmin
	pmat <- signif(-log10(pmat),3)
	ORmat <- log(enr$OR)
	heatmap.2(ORmat,Rowv=TRUE,Colv=TRUE,dendrogram="both", cellnote=pmat,
		notecol="black",col=redblue(256),scale="none",key=TRUE, keysize=1.5,
		density.info="none", trace="none", cexRow=0.5,cexCol=0.5,notecex=0.5)
}

## Defining a SubFunction for plotting a tSNE plot based on Module EigenGenes:
plot.ME.tSNE <- function(dat,net,lab,col){
	# Get Module EigenGenes:
	ME <- get.module.eg(dat,net)
	# Vizualize:
	require("Rtsne")
	tnse_out <- Rtsne(as.matrix(ME))
	plot(tsne_out2$Y,pch=22,bg=col,col="black")
}

