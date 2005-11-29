\name{sestck}
\alias{sestck}
\title{Tukey's biweight objectif function (scale)}
\description{
Computes Tukey's biweight objectif function (scale) corresponding 
with the Mahalanobis distances x.  
}
\usage{
sestck(x, start, c, k, tol)
}
\arguments{
  \item{x}{a vector of data corresponds to the Mahalanobis distances.}
  \item{start}{indicates the initial scale if different from zero, else scale will be estimated by the median.}
  \item{c}{Tukey's biweight constant.}
  \item{k}{the tuning constant for the multivariate S-estimator of location and covariance.}
  \item{tol}{indicates the tolerance value of the iteration algorithm.}
}
\keyword{multivariate}
\keyword{robust}

