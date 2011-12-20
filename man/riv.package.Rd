\name{riv-package}
\alias{riv-package}
\docType{package}
\title{Robust instrumental variables estimator package.
}
\description{
This package contains tools to find a robust instrumental variables estimator based on a high breakdown point S-estimator of location and
covariance.
}
\details{
\tabular{ll}{
Package: \tab riv\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2006-02-16\cr
License: \tab GLP\cr
}
slc(x, nsamp = 500, bdp = 0.5):

finds a high breakdown point S-estimator of multivariate location and covariance.

riv(Y,Xex=FALSE,Xend,W,intercept=TRUE,method=c("robust","classical"),nsamp=500,bdp=.5,approx=1000):

finds a robust instrumental variables estimator using a high breakdown point S-estimator of location and covariance.
}
\author{
Kaufmann B. \email{beat.kaufmann@epfl.ch}

Cohen-Freue G.V. \email{gcohen@stat.ubc.ca}

Zamar R.H. \email{ruben@stat.ubc.ca}
}
\references{
LOPUHAA,H.P. (1989). \emph{On the Relation between S-estimators and M-estimators of Multivariate Location and Covariance. Ann. Statist.} \bold{17} 1662-1683. 

COHEN-FREUE,G.V. and ZAMAR,R.H. (2005). \emph{A Robust Instrumental Variables Estimator.}
}
\keyword{multivariate}
\examples{
## simulation of a multivariate data-set
library(MASS)
x <- mvrnorm(100,c(0,0,0),matrix(c(1,0,0,0,1,0.5,0,0.5,1),ncol=3))

slc(x,nsamp=50,bdp=0.15)

## load data earthquake: the first column is the response Y, the second the endogenous variable X 
## and the third column is the instrument W.
data(earthquake)
riv(earthquake[,1],FALSE,earthquake[,2],earthquake[,3]) 
}
