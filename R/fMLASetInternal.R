#internal functions imported from LiquidAssociation package with permission
		stand <- function(object) {
			if (sum(is.na(object))>0)
			stop("Input vector contains missing value!!")
			myvect <- object
			ans <- (myvect-mean(myvect))/sd(myvect)
			return(ans)
		}
		qqnorm2 <- function(object){
			if (sum(is.na(object))>0)
			stop("Input vector contains missing value!!")
			myvect <- object
			n <- length(myvect)
			rmyvect <- rank(myvect)/(n+1)
			nmyvect <- sapply(rmyvect, qnorm)
			return(nmyvect)
		}
		permute3.se <- function(object){
			stand.data <- object
			null <- sample(stand.data[,3])	
			ndata <- cbind(stand.data[,1:2], null)
			return(ndata)
		}
		boots.gla.se2 <- function(object, num=1 ,cut=4, dim=3){
			stand.data <- object
			n <- nrow(stand.data)
			sam <- sample(1:n, n, replace=TRUE)
			ndata <- stand.data[sam,]
			ndata[,dim] <- qqnorm2(ndata[,dim])
			stand.ndata <- apply(ndata, 2, stand)
			ans <- GLA(stand.ndata, cut=cut, dim=dim)
			return(ans)
		}

#internal functions
exp.fun <- function(x){
	zero.one <- as.numeric(cut2(x,g=3))
	return(zero.one)
}

stand2 <- function(object){
	ans <- (object-mean(object,na.rm=TRUE))/sd(object,na.rm=TRUE)
	return(ans)
}

quant.norm <- function(vector){
	myvect <- vector
	n <- sum(as.numeric(!is.na(myvect)))
	rmyvect <- rank(myvect,na.last="keep")/(n+1)
	nmyvect <- sapply(rmyvect, qnorm)
	return(nmyvect)
}

# updates LiquidAssociation's GLA function to call Hmisc's cut2 vs base cut to avoid breaks are not unique error 
GLA <- function(object, cut=4, dim=3, geneMap=NULL){
	if (!is.numeric(object))
		stop("Input matrix must be numeric")
	if(ncol(object)!=3)
		stop("Input data must have three variables")
	if (length(colnames(object))!=3 & length(geneMap)!=3)
		stop("Please specify the names of three variables")
	if (length(geneMap)==3)
		colnames(object)<-names(geneMap)
	data <- object[!is.na(object[,1]) & !is.na(object[,2]) & !is.na(object[,3]),]
	data[,dim] <- qqnorm2(data[,dim])
	data <- apply(data, 2, stand)
	x <- seq(0,1, length=cut)
	br <- quantile(data[,dim], prob=x)
	index <- as.numeric(Hmisc::cut2(data[, dim], g=(cut-1)))
	tab <- table(index)
	tab <- tab[tab > 2]
	vect <- as.numeric(names(tab))
	m2 <- rep(0, length(vect))
	cor.e <- rep(0, length(vect))
	gla.x2 <- rep(0, length(vect))
	for ( i in 1:length(vect)){
		p <- which(index==vect[i])
		m2[i] <- mean(data[p,dim])
		cor.e[i] <- cor(data[p, -dim])[1,2]
		gla.x2[i] <- cor.e[i]*m2[i]
	}
	ans <- mean(gla.x2)
	names(ans) <- paste("GLA(", colnames(object)[1], ",", colnames(object)[2], "|", colnames(object)[3],")", sep="")
	return(ans)
	}

#####internal to fastMLA
jobsplit <- function(ival=1, data, nvec=NULL, rvalue=0.5, cut=4){
	#top <- matrix(NA,ncol=5,nrow=topn)
	top = data.frame(
		X       = character(),
		Y       = character(),
		Z       = character(),
		rhoDiff = double(),
		LA      = double()
	)
	#data.cor <- data[,-ival]
	third    <- data[,ival,drop = F]
	data.cor <- data[,-nvec$z,drop = F]
	# pwwang
	allnames = colnames(data)
	names3   = allnames[ival]
	if (is.null(nvec$x)) {
		names1 = allnames
		names2 = allnames
	} else {
		names1 = allnames[nvec$x]
		names2 = setdiff(allnames, c(names1, allnames[nvec$z]))
	}
	if(length(unique(third))<3 | nlevels(Hmisc::cut2(data[,ival], g=3))< min(3, cut-1)){
		return(top)
	} else {
		#data sorted based on column input, separates into high vs low expression 
		#based on gene in third position
		subset   <- !is.na(third)
		exp.vec  <- exp.fun(third)
		invec    <- exp.vec[subset]
		data.mat <- data.cor[subset,]
		ones     <- data.mat[which(invec==max(invec)),]
		zeroes   <- data.mat[which(invec==min(invec)),]
		#calculates correlation with call to WGCNA cor()
		#use="p" for pearson pairwise handling of missing data
		rho.one  <- WGCNA::cor(ones, use="p")
		rho.zero <- WGCNA::cor(zeroes, use="p")
		rho.diff <- rho.one-rho.zero
		rho.diff[upper.tri(rho.diff,diag=TRUE)] <- NA
		#determines which pairs are greater than the specified correlation
		#code to examine rhodiff values
		dnames = colnames(data.mat)
		GLAcalc <- function(subord,data.cor,third){
			trip <- cbind(data.cor[,subord],third)
			outGLA <- GLA(trip, dim=3, cut=cut)
			return(outGLA)
		}
		for(i in 1:(ncol(data.mat)-1)){
			# only keep res.mat[,1] as the X
			if (!dnames[i] %in% names1) next
			index  = which((abs(rho.diff[,i]))>rvalue, arr.ind=TRUE)
			inames = intersect(dnames[index], names2)
			index  = match(inames, dnames)
			if (length(inames)!=0){
				res.mat <- matrix(NA,ncol=5,nrow=length(index))
				colnames(res.mat) = colnames(top)
				res.mat[,1] <- dnames[i]	
				res.mat[,2] <- inames
				res.mat[,3] <- names3
				res.mat[,4] <- rho.diff[,i][index]
				if(length(inames)==1){
					res.mat[,5] = GLAcalc(res.mat[,1:2],data.cor=data.cor,third=third)
				} else {
					GLAtest <- apply(res.mat[,1:2],1,GLAcalc,data.cor=data.cor,third=third)
					res.mat[,5] = GLAtest
				}
				top <- rbind(top,res.mat)
			} else {top <- top}
		} 
		return(top)
	}
}
wrapper<-function(data, topn, nvec, rvalue, cut){
	outlist = mclapply(nvec$z, jobsplit, data=data, nvec = nvec, rvalue=rvalue, cut=cut)
	outlist = do.call(rbind,outlist)
	outlist = outlist[order(abs(as.numeric(outlist[,5])),decreasing=TRUE),,drop = F]
	if (!is.null(topn))
		outlist = head(outlist, topn)
	return(outlist)
}
namesfun<-function(toplt,data,nvec){
	#ival <- toplt[3]
	long  <- colnames(data)
	if (is.null(nvec$x)) {
		shortx <- long[-nvec$z]
		shorty <- NULL
	} else {
		shortx <- long[nvec$x]
		shorty <- setdiff(long, shortx, long[topIt[3]])
	}
	outname <- cbind(shortx[toplt[1]],ifelse(is.null(shorty), shortx, shorty)[toplt[2]],long[toplt[3]])
	return(outname)
}
###end of internal to fastMLA

####internal to mass.CNM
#for running samples with CNM.full in bulk to return top values
top.CNM <- function(data, GLA.mat){
	if (!is.numeric(data))
		stop("Data matrix must be numeric")
	dat.q <- apply(data,2,quant.norm)
	dat.s <- apply(dat.q,2,stand2)
	pvalmat <- matrix(NA,nrow=nrow(GLA.mat),ncol=10)
	pvalmat[,10] <- "F"
	pvalmat[,1:5] <- as.matrix(GLA.mat)
	for (i in 1:nrow(pvalmat)){
		datCNM <- dat.s[,pvalmat[i,1:3]]
		pvalmat[i,6:9]<-CNM.full(datCNM)@output[8,]
	}
	colnames(pvalmat) <- c("X1 or X2","X2 or X1","X3","rhodiff","GLA value","estimates","san.se","wald","p value","model")
	return(pvalmat)
}

#for running samples with CNM.simple in bulk
sens.CNM<-function(data, full.mat){
	if (!is.numeric(data))
		stop("Data matrix must be numeric")
	dat.q <- apply(data,2,quant.norm)
	dat.s <- apply(dat.q,2,stand2)
	pvalmat <- matrix(NA,nrow=nrow(full.mat),ncol=10)
	pvalmat[,10] <- "S"
	pvalmat[,1:5] <- as.matrix(full.mat[,1:5])
	for (i in 1:nrow(pvalmat)){
		datCNM <- dat.s[,pvalmat[i,1:3]]
		pvalmat[i,6:9] <- CNM.simple(datCNM)@output[4,]
	}
	colnames(pvalmat) <- c("X1 or X2","X2 or X1","X3","rhodiff","GLA value","estimates","san.se","wald","p value","model")
	return(pvalmat)
}

#used to determine if results from CNM.full or CNM.simple were sensible
#intended for use with output from sens.CNM or top.CNM
#to determine results which need to be run in next stage of pvalue testing
boots.index <- function(matrix){
	index <- ((as.numeric(matrix[,7])==0)|(as.numeric(matrix[,7])>10)|(is.finite(as.numeric(matrix[,6]))!=TRUE)|
		(is.finite(as.numeric(matrix[,7]))!=TRUE)|(is.finite(as.numeric(matrix[,8]))!=TRUE)|
		(is.finite(as.numeric(matrix[,9]))!=TRUE))
	return(index)
	}

####end of internal to mass.CNM

####internal to fastboots.GLA
##large cluster code for bootstrapping which avoids the chance NA iteration error
##cluster version rewrite of getGLA
clusterGLA <- function(object, boots=30, clust, perm=100, cut=4, dim=3, geneMap=NULL){
		if (!is.numeric(object))
			stop("Input matrix must be numeric")
		if (ncol(object)!=3)
			stop("Input data must have three variables")
		 if (length(colnames(object))!=3 & length(geneMap)!=3)
			stop("Please specify the names of three variables")
		if (length(geneMap)==3)
			colnames(object)<-names(geneMap)
		dat <- object[(!is.na(object[,1]) & !is.na(object[,2]) & !is.na(object[,3])),]
	#make X3 into its normal quantile and for all 3 transform to mean 0,sd 1
		dat[,dim] <- qqnorm2(dat[,dim])
		nsdata <- apply(dat, 2, stand)
	#obtain initial GLA estimate from transformed data
		gla1 <- GLA(nsdata, cut=cut, dim=dim)
	#boots.gla.se2 randomly samples from the transformed data, then re-quantile normalizes
	#only the X3 portion before obtaining a GLA estimate
	#using the sapply function, it obtains #=boots estimates of GLA in this way
		bootsse1 <- sapply(1:boots, boots.gla.se2, object=nsdata, cut=cut, dim=dim)
 		bse1 <- sd(bootsse1)
	#bse1 is sd of bootstrap results, test1 is statistic created by dividing original
	#GLA estimate by bootstrap obtained sd
		test1 <- gla1/bse1
	#nullGLA
	nullGLA <- function(iter,boots,data,cut,dim,clust=clust){
		#.libPaths(.libPaths())
		#library('LiquidAssociation')
		stand <- function(object) {
			if (sum(is.na(object))>0)
			stop("Input vector contains missing value!!")
			myvect <- object
			ans <- (myvect-mean(myvect))/sd(myvect)
			return(ans)
		}
		qqnorm2 <- function(object){
			if (sum(is.na(object))>0)
			stop("Input vector contains missing value!!")
			myvect <- object
			n <- length(myvect)
			rmyvect <- rank(myvect)/(n+1)
			nmyvect <- sapply(rmyvect, qnorm)
			return(nmyvect)
		}
		permute3.se <- function(object){
			stand.data <- object
			null <- sample(stand.data[,3])	
			ndata <- cbind(stand.data[,1:2], null)
			return(ndata)
		}
		boots.gla.se2 <- function(object, num=1 ,cut=4, dim=3){
			stand.data <- object
			n <- nrow(stand.data)
			sam <- sample(1:n, n, replace=TRUE)
			ndata <- stand.data[sam,]
			ndata[,dim] <- qqnorm2(ndata[,dim])
			stand.ndata <- apply(ndata, 2, stand)
			ans <- GLA(stand.ndata, cut=cut, dim=dim)
			return(ans)
		}
		GLA <- function(object, cut=4, dim=3, geneMap=NULL){
		if (!is.numeric(object))
			stop("Input matrix must be numeric")
		if(ncol(object)!=3)
			stop("Input data must have three variables")
		if (length(colnames(object))!=3 & length(geneMap)!=3)
			stop("Please specify the names of three variables")
		if (length(geneMap)==3)
			colnames(object)<-names(geneMap)
		data <- object[!is.na(object[,1]) & !is.na(object[,2]) & !is.na(object[,3]),]
		data[,dim] <- qqnorm2(data[,dim])
		data <- apply(data, 2, stand)
		x <- seq(0,1, length=cut)
		br <- quantile(data[,dim], prob=x)
		index <- as.numeric(Hmisc::cut2(data[, dim], g=(cut-1)))
		tab <- table(index)
		tab <- tab[tab > 2]
		vect <- as.numeric(names(tab))
		m2 <- rep(0, length(vect))
		cor.e <- rep(0, length(vect))
		gla.x2 <- rep(0, length(vect))
		for ( i in 1:length(vect)){
			p <- which(index==vect[i])
			m2[i] <- mean(data[p,dim])
			cor.e[i] <- cor(data[p, -dim])[1,2]
			gla.x2[i] <- cor.e[i]*m2[i]
		}
		ans <- mean(gla.x2)
		names(ans) <- paste("GLA(", colnames(object)[1], ",", colnames(object)[2], "|", colnames(object)[3],")", sep="")
		return(ans)
	 	}
		datanull <- permute3.se(data)
		gla2 <- GLA(datanull,cut=cut,dim=dim)
		bootsse2 <- sapply(1:boots,boots.gla.se2,object=datanull,cut=cut,dim=dim)
		bse2 <- sd(bootsse2)
		ans <- gla2/bse2
		return(ans)}
		iter <- perm
		check <- parallel::parLapply(clust,1:iter,nullGLA,boots=boots,data=nsdata,cut=cut,dim=dim)
		ans <- unlist(check)
	#the pvalue is 2*# of times permuted statistic exceeds abs(original GLA)
	#framing the p-value in this way avoids the errors caused by NA through sampling
		pvalue <- 2*sum(abs(test1) <= ans,na.rm=TRUE)/sum(as.numeric(!is.na(ans)))
		out <- c(test1, pvalue)
		names(out) <- c("sGLA", "p value")
		return(out)
		}

#creates list of triplets from data based on list of results to bootstrap
makelist <- function(matrix,data){
testlist <- list()
for(i in 1:nrow(matrix)){
testlist[[i]] <- data[,matrix[i,1:3]]
}
return(testlist)
}
#####end of internal to fastboots.GLA
