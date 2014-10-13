\name{fastboots.GLA-methods}
\docType{methods}
\alias{fastboots.GLA-methods}
\alias{fastboots.GLA-class}
\alias{fastboots.GLA,ANY,matrix-method}
\alias{fastboots.GLA}
\title{Function to parallelize bootstrapping for MLA robust estimates.}
\description{
Function to increase speed of bootstrapping for MLA robust estimates. It performs the bootstrapping in parallel, so the time decrease will be dependent on the number of cores or CPUs available. Intended for use with results of triplets not processed sensibly by call to either CNM.full or CNM.simple in mass.CNM, though it can also be used with the direct results of fastMLA.
}
\usage{
fastboots.GLA(tripmat, data, clust, boots = 30, perm = 100, cut = 4)
}
\arguments{
  \item{tripmat}{
Matrix specifying the triplets to be estimated. Intended for use with the results of the mass.CNM or fastMLA functions.
}
  \item{data}{
Matrix of numeric data, with columns representing genes and rows representing observations.
}
  \item{clust}{
Cluster of CPU cores created by makeCluster. See parallel.
}
  \item{boots}{
Number of bootstrap iterations to estimate direct estimate SE.
}
  \item{perm}{
Number of permutations to use in estimating p-value.
}
  \item{cut}{
Value passed to the internal clusterGLA function to create buckets (equal to number of buckets+1). Values placing between 15-30 samples per bucket are optimal. Must be a positive integer>1. See GLA.
}
}
\details{
Choosing the number of bins:
For example, assume that our data has 100 observations. Since values between 15-30 observations per bin are optimal, good values to choose would be 5-7.
}
\value{
fastboots.GLA returns a data.frame that specifies the genes in positions X1, X2 and X3; the rhodiff value of the triplet, the GLA value of the triplet, the direct estimate statistic, and the direct estimate p-value. More detailed explanation is available in the package vignette. 
}
\keyword{methods}
\references{
[1] Yen-Yi Ho, Giovanni Parmigiani, Thomas A Louis, and Leslie M Cope. Modeling liquid association.
Biometrics, 67(1):133-141, 2011.
}
\author{
Tina Gunderson
}

\section{Warning }{
The data matrix must be numeric. 
} 

\note{
The cluster of CPUs to use for bootstrapping must be created with makeCluster from the parallel package before running the fastboots.GLA function. 
}

\seealso{
fastMLA, mass.CNM, LiquidAssociation, parallel
}
\examples{
#to view the code for the function
selectMethod("fastboots.GLA", signature=c("ANY","matrix"))

#
library(yeastCC)
data(spYCCES)
lae <- spYCCES[,-(1:4)]
### get rid of samples with high NA elements
lae <- lae[apply(is.na(exprs(lae)),1,sum) < ncol(lae)*0.3,]
data <- t(exprs(lae))
data <- data[,1:50]
dim(data)

example <- fastMLA(data=data, topn=25, nvec=1:10, rvalue=1.0, cut=4)
clust <- makeCluster(4)
ex <- example[1:5,]
GLAeasy <- fastboots.GLA(ex, data=data, clust=clust, boots=30, perm=100, cut=4)
GLAeasy
stopCluster(clust)
closeAllConnections()
}
