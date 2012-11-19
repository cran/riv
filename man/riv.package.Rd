\name{riv-package}
\alias{riv-package}
\docType{package}
\title{Robust instrumental variables estimator package.}
\description{
  This package contains tools to find a robust instrumental variables
  estimator based on a high breakdown point S-estimator of location and
  covariance.
}
\details{
  \tabular{ll}{
    Package: \tab riv\cr
    Type: \tab Package\cr
    Version: \tab 2.0\cr
    Date: \tab 2012-11-19\cr
    License: \tab GPL-2\cr
  }
  riv(Y, Xend, Xex = NULL, Zinst, dummies = NULL,
      method = c('S-est', 'SD-est', 'MCD-est', 'classical'),
      nsamp = 500, bdp = 0.5)

  finds a robust instrumental variables estimator using a high breakdown
  point S-estimator of location and covariance.
}
\author{
  G.V. Cohen-Freue \email{gcohen@stat.ubc.ca}

  D. Cubranic \email{cubranic@stat.ubc.ca}
  
  with contributions from B. Kaufmann \email{beat.kaufmann@epfl.ch} and
  R.H. Zamar \email{ruben@stat.ubc.ca}
}
\references{
LOPUHAA,H.P. (1989). \emph{On the Relation between S-estimators and M-estimators of Multivariate Location and Covariance. Ann. Statist.} \bold{17} 1662-1683. 

COHEN-FREUE, G.V., ORTIZ-MOLINA, H., and ZAMAR, R.H. (2012) \emph{A
  Natural Robustification of the Ordinary Instrumental Variables
  Estimator.} Submitted to \emph{Biometrics}.
}

\keyword{multivariate}

\examples{
## load data earthquake: the first column is the response Y, the second
## the endogenous variable X and the third column is the instrument
## Zinst.
data(earthquake)
riv(earthquake[,1], earthquake[,2], NULL, earthquake[,3])
}
