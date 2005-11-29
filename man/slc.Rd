\name{slc}
\alias{slc}
\title{S-estimator of multivariate location and covariance}
\description{Finds the robust S-estimator of multivariate location and covariance with a
high breakdown point based on Tukey's biweight function.
}
\usage{
slc(x, nsamp = 500, bdp = 0.5)
}

\arguments{
  \item{x}{a matrix or data frame of data values, say of dimension n x p. The rows represent observations and the columns represent variables.}
  \item{nsamp}{the number of random p-subsets considered to compute the S-estimator (default = 500).}
  \item{bdp}{breakdown point value of the S-estimator must be 0.15,0.25 od 0.5 (default).}
}
\value{
  A list with components
  \item{location}{vector of the estimated multivariate location.}
  \item{covariance}{matrix of the estimated covariance.}
  \item{distances}{vector of robust Mahalanobis distances versus location and covariance.}
  \item{scale}{distance of the scale estimates.}
  \item{c}{constant of the Tukey's biweight function.}
}
\references{
LOPUHAÄ,H.P. (1989) \emph{On the Relation between S-estimators and M-estimators of Multivariate Location and Covariance. Ann. Statist.} \bold{17} 1662-1683.

RUPPERT,D. (1992). \emph{Computing S Estimators for Regression and Multivariate Location/Dispersion. J. Comput. Graph. Statist.} \bold{1} 253-270.
}
\author{
Kaufmann B. \email{beat.kaufmann@epfl.ch} 

Cohen-Freue G.V. \email{gcohen@stat.ubs.ca}
}

\seealso{\code{\link{lqs}}}

\examples{ 
library(MASS)
## load the mortality data-set (62 Alaskan earthquake observations)
data(mortality)
slc.res <- slc(mortality,nsamp=100,bdp=0.25)   

## plot of the Mahalanobis Distances
plot(slc.res[[3]],xlab="Observation",ylab="Mahalanobis Distance",main="The Mahalanobis Distance of each observation")   

## simulation of a multivariate data-set
library(MASS)
x <- mvrnorm(100,c(0,0,0),matrix(c(1,0,0,0,1,0.5,0,0.5,1),ncol=3))

slc(x,nsamp=50,bdp=0.15)
}
\keyword{robust}
\keyword{multivariate}
