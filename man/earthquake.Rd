\name{earthquake}
\alias{earthquake}
\docType{data}
\title{Alaskan earthquakes data}
\description{
The dataset contains information about 62 Alaskan earthquakes occured from the period of 1989-1978.
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
FULLER,W.A. (1987). \emph{Measurement Error Models}. Wiley, NewYork. 
}
\references{
COHEN-FREUE,G.V. and ZAMAR,R.H. (2005). \emph{A Robust Instrumental Variables Estimator}.
}
\examples{
data(earthquake)
plot(earthquake[,2],earthquake[,1],xlab="X",ylab="Y") 
}
\keyword{datasets}
