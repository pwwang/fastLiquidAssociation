\name{mass.CNM-methods}
\docType{methods}
\alias{mass.CNM-methods}
\alias{mass.CNM-class}
\alias{mass.CNM,matrix,data.frame-method}
\alias{mass.CNM}
\title{Function to obtain CNM significance estimates from an object result of fastMLA.}
\description{
Function which takes the results of fastMLA and estimates significance based on the beta5 values returned from the conditional normal model estimate. As there are two conditional normal model functions available, full and simple (see CNM.simple, CNM.full, or the fastLiquidAssociation vignette for further details), it attempts to first process all triplets through the full model and for those that do not appear to be well fit by the full model, processes those triplets through the simple model. 
It returns the triplet(s), fastMLA correlation difference, MLA estimate, beta5 results from the CNM.full or CNM.simple function, and model type, sorted by most significant results according beta5 p-value. Also returns a list of triplets not sensibly fit by either model.
}
\usage{
mass.CNM(data, GLA.mat, nback = 100)
}
\arguments{
  \item{data}{
Matrix of numeric data, with columns representing genes and rows representing observations.
}
  \item{GLA.mat}{
Matrix of triplets to be processed. Intended for use with results of fastMLA function. The 	X3 position will be defaulted to the value in the third column. see fastMLA.
}
  \item{nback}{
Number of results to return. Results are sorted from most significant p-value. 
}
}
\value{
\item{top p-values}{A data frame with 10 variables: the genes in positions X1, X2 and X3; the rhodiff value of the triplet; the GLA value of the triplet; the estimates for the beta_5 value, standard error, Wald test statistic, and p-value for the Wald test statistic; and the conditional normal model the estimates were obtained from (either F or S, standing for full and simple respectively). A more detailed discussion of these values is available in the vignette.}
\item{bootstrap triplets}{A matrix with triplets that did not appear to be well fit by either the full or simple CNM model and thus are recommended to have significance estimated via fastboots.GLA. Its variables are otherwise as specified in the top p-values results description.}
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
Observing the warning message: "In sqrt(diag(object$valpha)) : NaNs produced" is not uncommon in the mass.CNM function as it can be caused by failures to fit either the full or simple CNM model to a triplet.
} 

\seealso{
fastMLA, LiquidAssociation
}
\examples{
#to view function code
selectMethod("mass.CNM", signature=c("matrix","data.frame"))

#
library(LiquidAssociation)
library(yeastCC)
data(spYCCES)
lae <- spYCCES[,-(1:4)]
### get rid of samples with high % NA elements
lae <- lae[apply(is.na(exprs(lae)),1,sum) < ncol(lae)*0.3,]
data <- t(exprs(lae))
data <- data[,1:50]
example <- fastMLA(data=data, topn=25, nvec=1:10, rvalue=1.0, cut=4)
CNMcalc <- mass.CNM(data=data, GLA.mat=example, nback=10)
CNMcalc
}
