\name{earthquake}
\alias{earthquake}
\docType{data}
\title{Alaskan earthquakes data}
\description{

  The dataset contains information about 62 Alaskan earthquakes that
  occured between 1969-1978 (Fuller, 1987). The goal is to see how the
  earthquake strength, measured in terms of the true value of the body
  waves, \eqn{x^{*}}, impacts on the amplitude of the surface waves of
  the earthquake (\eqn{Y}). However, we do not observe \eqn{x^{*}} but
  \eqn{X}, which is the logarithm of the seismogram amplitude of
  longitudinal body waves measured at some observation stations, i.e.,
  \eqn{x_i=x^{*}_i+u_i}, for \eqn{i=1,\dots,62}, where \eqn{u_i} is the
  measurement error of each observation. Thus, in the regression
  \eqn{Y_i=a+b*X_i+e_i}, the covariate \eqn{X} is an endogenous
  covariate. We can consistently estimate the regression parameters
  using instrumental variables estimators. The logarithm of maximum
  seismogram trace amplitude at short distance, \eqn{W}, can be used as
  an instrument.

  The first column in the dataset is the response (Y), the second column
  is the endogenous variable (X), and the third column is the instrument
  (W).
}

\usage{data(earthquake)}
\format{
  A data frame with 62 observations on the following 3 variables.
  \describe{
    \item{\code{Y}}{a numeric vector of the logarithm of the seismogram amplitude of 20 second waves.}
    \item{\code{X}}{a numeric vector of the logarithm of the seismogram amplitude of longitudinal body waves.}
    \item{\code{W}}{a numeric vector of the logarithm of maximum seismogram trace amplitude at short distance.}
  }
}
\source{
FULLER,W.A. (1987). \emph{Measurement Error Models}. Wiley, New York. 
}
\references{
COHEN-FREUE,G.V. and ZAMAR,R.H. (2005). \emph{A Robust Instrumental Variables Estimator}.
}
\examples{
data(earthquake)
plot(earthquake$X,earthquake$Y,xlab="X",ylab="Y") 
}
\keyword{datasets}
