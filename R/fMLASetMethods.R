#################
###Method setting

#####for fastMLA
fastMLA <- setClass("fastMLA",slots=c(fastMLA="data.frame"))

setGeneric("fastMLA", function(data,topn=NULL,nvec=1, rvalue=0.5,cut=4,zcat=F,threads=detectCores()) standardGeneric("fastMLA"))

setMethod("fastMLA",signature(data="matrix"),
	function(data,topn=NULL,nvec=1, rvalue=0.5,cut=4,zcat=F,threads=detectCores()){
		if(ncol(data)<3)
			stop("There must be at least 3 columns in the data set")
		if(!is.list(nvec)) nvec = list(z = nvec)
		if(ncol(data)<length(nvec$z))
			stop("Specified number of genes to test exceeds number of genes in data matrix") 
		if(!is.numeric(data))
			stop("Data matrix must be numeric")
		if(!is.null(topn)&&(topn<1||!is.numeric(topn)||(topn%%1!=0)))
			stop("topn must be a positive whole number")
		if(any(nvec$z<1)|any(!is.numeric(nvec$z))|any(nvec$z%%1!=0))
			stop("nvec must contain positive whole number(s)")
		if(rvalue<0|rvalue>2)
			stop("rvalue must be between 0 and 2")
		if(cut<0|!is.numeric(cut))
			stop("cut must be a positive whole number")
		if (threads > 1) enableWGCNAThreads(threads)
		dat.n = data
		if (zcat) dat.n = data[, -nvec$z, drop=F]
		dat.q <- apply(dat.n,2,quant.norm)
		dat.s <- apply(dat.q,2,stand2)
		if (zcat) {
			data[, -nvec$z] = dat.s
		} else {
			data = dat.s
		}
		return(wrapper(data,topn,nvec,rvalue,cut,zcat))
		if (threads > 1) doParallel::stopImplicitCluster()
	}
)

#####for mass.CNM
mass.CNM <- setClass("mass.CNM",slots=c(mass.CNM="list"))

setGeneric("mass.CNM", function(data,GLA.mat, nback=100) standardGeneric("mass.CNM"))

setMethod("mass.CNM",signature(data="matrix",GLA.mat="data.frame"),
	function(data, GLA.mat, nback=100){
	if(ncol(data)<3)
		stop("There must be at least 3 columns in the data set")
	if(!is.numeric(data))
		stop("Data matrix must be numeric")
	if(ncol(GLA.mat)<5)
		stop("GLA.mat must have 5 columns")
	if(nback<1|!is.numeric(nback)|any(nback%%1!=0))
		stop("nback must be a positive integer")
	if(any(nback<2))
		stop("nback must be >= 2")
	fullmod <- top.CNM(data, GLA.mat)
	boots.f <- boots.index(fullmod)
	rpts.full <- fullmod[boots.f,]
	comp.full <- fullmod[which(boots.f==FALSE),]
	rpts.simple <- NULL
	if(sum(as.numeric(boots.f))!=0){
	rpts.full<-matrix(rpts.full,ncol=10)
	simpmod <- sens.CNM(data=data, full.mat=rpts.full)
	boots.s <- boots.index(simpmod)
	rpts.simple <- simpmod[boots.s,]
	rpts.simple[,4:9] <- round(apply(rpts.simple[,4:9],2,as.numeric),4)
	comp.simp <- simpmod[which(boots.s==FALSE),]
	pvalmat <- rbind(comp.full,comp.simp)
	highps <- head(pvalmat[order(as.numeric(pvalmat[,9]),-(as.numeric(pvalmat[,8])),decreasing=F),],min(nback,nrow(pvalmat)))
	highps[,4:9] <- round(apply(highps[,4:9],2,as.numeric),4)
	} else {
	pvalmat <- comp.full
	highps <- head(pvalmat[order(as.numeric(pvalmat[,9]),-(as.numeric(pvalmat[,8])),decreasing=FALSE),],min(nback,nrow(pvalmat)))
	highps[,4:9] <- round(apply(highps[,4:9],2,as.numeric),2)
	}
	convert <- matrix(as.numeric(highps[,4:9]),byrow=FALSE,ncol=6)
	highps <- data.frame(highps[,1:3],convert,highps[,10],stringsAsFactors=FALSE)
	colnames(highps) <- c("X1 or X2","X2 or X1","X3","rhodiff","MLA value","estimates","san.se","wald","p value","model")
	output <- list(highps, rpts.simple)
	names(output) <- c("top p-values","bootstrap triplets")
	return(output)
	stopImplicitCluster()
	}
)


#####for fastboots.GLA
fastboots.GLA <- setClass("fastboots.GLA",slots=c(fastboots.GLA="data.frame"))

setGeneric("fastboots.GLA", function(tripmat,data,clust,boots=30, perm=100, cut=4) standardGeneric("fastboots.GLA"))

###setMethod("fastboots.GLA", signature(tripmat="list"), 
###	function(tripmat,data,...){
###		tripmat<-tripmat[[1]]
###		fastboots.GLA(tripmat,data,...)
###	}
###)

setMethod("fastboots.GLA",signature(data="matrix"),
	function(tripmat,data,clust,boots=30, perm=100, cut=4){
	if(ncol(data)<3)
		stop("There must be at least 3 columns in the data set")
	if(!is.numeric(data))
		stop("Data matrix must be numeric")
	if(!is.data.frame(tripmat)&!is.matrix(tripmat))
		stop("tripmat must be a data frame or matrix")
	if(ncol(tripmat)<5)
		stop("tripmat must have 5 columns")
	if(any(boots<1)|any(!is.numeric(boots))|any(boots%%1!=0))
		stop("boots must contain positive whole number(s)")
	if(any(perm<1)|any(!is.numeric(perm))|any(perm%%1!=0))
		stop("perm must contain positive whole number(s)")
	if(any(cut<2)|any(!is.numeric(cut))|any(cut%%1!=0))
		stop("cut must contain positive whole number(s)")
	dat.q <- apply(data,2,quant.norm)
	dat.s <- apply(dat.q,2,stand2)
	bootslist <- NULL
	if (is.numeric(nrow(tripmat))==FALSE) {
	bootslist[[1]] <- dat.s[,tripmat[1:3]]
	tripmat <- matrix(tripmat,nrow=1)
	resmat <- clusterGLA(bootslist[[1]],boots=boots, clust, perm=perm, cut=cut, dim=3)
	results <- matrix(unlist(resmat),ncol=2,byrow=TRUE)
	} else {
	tripmat <- as.matrix(tripmat)
	bootslist <- makelist(tripmat,dat.s)
	resmat <- lapply(bootslist,clusterGLA,boots=boots, clust, perm=perm, cut=cut, dim=3)
	results <- matrix(unlist(resmat),ncol=2,byrow=TRUE)
	}
	newmat <- data.frame(matrix(NA,nrow=nrow(results),ncol=7))
	newmat[,1:5] <- tripmat[,1:5]
	newmat[,6:7] <- results
	colnames(newmat) <- c("X1 or X2","X2 or X1","X3","rhodiff","MLA value","MLA stat","boots p-value")
	#stopCluster(clust)
	return(newmat)
	stopImplicitCluster()
	}
)
