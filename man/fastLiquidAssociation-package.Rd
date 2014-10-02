\name{fastLiquidAssociation-package}
\alias{fastLiquidAssociation-package}
\alias{fastLiquidAssociation}
\docType{package}
\title{
Functions to extend liquid association (LA) analysis for genome-wide use
}
\description{
This package provides three external functions to extend liquid association (LA) analysis for genome-wide use: fastMLA, mass.CNM, and fastboots.GLA. It contains updated versions of several functions available in the LiquidAssociation package, but whereas those functions are written to be used with a single triplet, this package's functions are intended to be applied against a full data set.
}
\details{
\tabular{ll}{
Package: \tab fastLA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-08\cr
License: \tab GPL-2\cr
}
Contains three external functions: fastMLA, mass.CNM, and fastboots.GLA. fastMLA uses an algorithm to approximate liquid association values in order to screen for a subset likely to have high LA values and then evaluates that subset for their LA values. mass.CNM takes the results from fastMLA and attempts to estimate significance based on conditional normal models. In the event that mass.CNM is unable to provide an estimate, fastboots.GLA parallelizes the bootstrapping required to produce a more robust direct estimate for LA. 
}
\author{
Author: Tina Gunderson
Maintainer: Tina Gunderson <gunderson.tina@gmail.com>
}
\references{
[1] Yen-Yi Ho, Giovanni Parmigiani, Thomas A Louis, and Leslie M Cope. Modeling liquid association.
Biometrics, 67(1):133-141, 2011.

[2] Ker-Chau Li. Genome-wide coexpression dynamics: theory and application. Proceedings of the National
Academy of Sciences, 99(26):16875-16880, 2002.

[3] Paul T Spellman, Gavin Sherlock, Michael Q Zhang, Vishwanath R Iyer, Kirk Anders, Michael B Eisen,
Patrick O Brown, David Botstein, and Bruce Futcher. Comprehensive identification of cell cycle regulated
genes of the yeast saccharomyces cerevisiae by microarray hybridization. Molecular biology of the cell,
9(12):3273-3297, 1998.
}

\keyword{ package }
\seealso{
LiquidAssociation,WGCNA,parallel,yeastCC,GOstats,org.Sc.sgd.db
}
\examples{
library(fastLiquidAssociation)
library(yeastCC)
data(spYCCES)
lae <- spYCCES[,-(1:4)]
### get rid of samples with high % NA elements
lae <- lae[apply(is.na(exprs(lae)),1,sum) < ncol(lae)*0.3,]
data <- t(exprs(lae))
data <- data[,1:50]

##fastMLA
example <- fastMLA(data=data, topn=25, nvec=1:10, rvalue=1.0, cut=4, threads = detectCores())
example[1:5,]

##mass.CNM
CNMcalc <- mass.CNM(data=data, GLA.mat=example, nback=10)
CNMcalc

##fastboots.GLA
clust <- makeCluster(4)
ex <- example[1:5,]
GLAeasy <- fastboots.GLA(ex, data=data, clust=clust, boots=30, perm=100, cut=4)
GLAeasy
stopCluster(clust)
closeAllConnections()
}
 