\name{mortality}
\alias{mortality}
\docType{data}
\title{Mortality rate data}
\description{
The dataset contains information about the mortality rate from 60 U.S. cities using aggregate information from the year 1969/70. 
}
\usage{data(mortality)}
\format{
  A data frame with 60 observations on the following 8 variables.
  \describe{
    \item{\code{MO70}}{a numeric vector of the total mortality (number of deaths per 1000 people) from 1970.}
    \item{\code{MAGE}}{a numeric vector of the median age of the population (in years) from 1969.}
    \item{\code{CI68}}{a numeric vector of the number of packs of cigarettes per year per person.}
    \item{\code{MDOC}}{a numeric vector of the density of medical doctors (number of medical doctors per 100,000 people).}
    \item{\code{DENS}}{a numeric vector of the percentage of households with more than 1.5 persons per room.}
    \item{\code{NONW}}{a numeric vector of the fraction of the non-white population.}
    \item{\code{EDUC}}{a numeric vector of the percentage of the population over age 25 having a high-school diploma.}
    \item{\code{IN69}}{a numeric vector of the median income from 1969.}
  }
}

\source{
 CROCKER,D.T. et al. (1979). \emph{Methods
  Development for Assessing Air Pollution Control Benefits, Vol. 1.
  Experiments in the Economics of Epidemiology.} EPA-600/5-79-001a.
  Springfield, VA; National Technical Information Service. 
}
\examples{
data(mortality)
}
\keyword{datasets}
