\name{fastMLA-methods}
\docType{methods}
\alias{fastMLA-methods}
\alias{fastMLA-class}
\alias{fastMLA,matrix-method}
\alias{fastMLA}
\title{Function to more efficiently screen for gene triplets for those with a high liquid association value.}
\description{
Function reduces the processing power and memory needed to calculate modified liquid association (MLA) values for a genome by using a pre-screening method to reduce the candidate pool to triplets likely to have a high MLA value. It does this using matrix algebra to create an approximation to the direct MLA estimate for all possible pairs of X1X2|X3. 
}
\usage{
fastMLA(data, topn = 2000, nvec = 1, rvalue = 0.5, cut = 4, threads = detectCores())
}
\arguments{
  \item{data}{
Matrix of numeric data, with columns representing genes and rows representing observations.
}
  \item{topn}{
Number of results to return, ordered from highest |MLA| value descending. 
}
  \item{nvec}{
Numeric vector of the gene(s) to use in the X3 position of the X1X2|X3 screening. This should be a numeric vector representing the column #(s) of the gene.
}
  \item{rvalue}{
Tolerance value for LA approximation. Lower values of rvalue will cause a more thorough search, but take longer.
}
  \item{cut}{
Value passed to the GLA function to create buckets (equal to number of buckets+1). Values placing between 15-30 samples per bucket are optimal. Must be a positive integer>1. See GLA.
}
  \item{threads}{
Number of cores to use for multi-threading in correlation calculation (enableWGCNAThreads argument). See WGCNA.
}
}

\details{
Choosing the number of bins:
For example, assume that our data has 100 observations. Since values between 15-30 observations per bin are optimal, good values to choose for cut would be 5-7.
}
\value{
A data frame with 5 variables: the genes in positions X1, X2 and X3; the rhodiff value of the triplet; and the GLA value of the triplet. A more comprehensive discussion of these values is available in the vignette.
}
\references{
[1] Yen-Yi Ho, Giovanni Parmigiani, Thomas A Louis, and Leslie M Cope. Modeling liquid association.
Biometrics, 67(1):133-141, 2011.
}
\author{
Tina Gunderson
}

\note{
While this is intended to significantly reduce processing time for identifying high MLA values (and in our estimates did so by >90% at the default rvalue), even in a smaller data set (N=6000 genes), a 90% reduction still requires ~1.8 million calculations per gene and is still relatively time and memory intensive. 
}

\section{Warning }{
The data matrix must be numeric. 
} 

\seealso{
\code{\link{LiquidAssociation}}, \code{\link{parallel}}, \code{\link{WGCNA}}
}

\examples{
#to view function code
selectMethod("fastMLA", "matrix")

#
library(fastLiquidAssociation)
library(yeastCC)
data(spYCCES)
lae <- spYCCES[,-(1:4)]
### get rid of samples with high % NA elements
lae <- lae[apply(is.na(exprs(lae)),1,sum) < ncol(lae)*0.3,]
data <- t(exprs(lae))
data <- data[,1:50]

example <- fastMLA(data=data, topn=25, nvec=1:10, rvalue=1.0, cut=4)
example[1:5,]
closeAllConnections()
}