\name{riv}
\alias{riv}
\title{Robust instrumental variables estimator}

\description{
  Finds robust instrumental variables estimator using high breakdown
  point S-estimator of multivariate location and covariance.
}

\usage{
riv(Y, Xend, Xex=NULL, Zinst, dummies=NULL,
    method = c('S-est', 'SD-est', 'MCD-est', 'classical'))
}

\arguments{
  \item{Y}{vector of responses.}
  \item{Xend}{matrix of the endogenous variables, i.e. covariates that
    are correlated with the error.
  }
  \item{Xex}{matrix of the exogenous variables, i.e. covariates that are
    uncorrelated with the error. Default = \code{NULL}, i.e. no
    exogenous variables are present in the model.
  }
  \item{Zinst}{matrix of instruments, variables correlated with the
    endogenous covariates, but uncorrelated with the error term. The
    number of instrumental variables needs to be larger than or equal to
    the number of endogenous covariates.}
  \item{dummies}{matrix of exogenous dummy covariates, i.e.,
    where each \eqn{D_i} are 0--1 vectors.
  }
  \item{method}{the method to be used. The "\code{S-est}" method
    (default) is based on the S-estimator of multivariate location and
    covariance, and "\code{classical}" method is based on the sample
    mean and covariance. See \code{Details} section.
  }
}

\details{
  For method "\code{S-est}", the estimator is constructed by using a
  robust S-estimator of multivariate location and covariance based on
  the Tukey's biweight function (see \code{CovSest}).
  
  If RIV is computed using the S-estimator, its variance-covariance
  matrix is estimated based on the empirical influence function. See
  references for more details.

  For method "\code{SD-est}", the estimator is constructed by using the
  Stahel-Donoho robust multivariate location and scatter estimator (see
  \code{CovSde}).
  
  For method "\code{MCD-est}", the estimator is constructed by using a
  'Fast MCD' (Minimum Covariance Determinant) robust multivariate
  location and scatter estimator (see \code{CovMcd}).
  
  For method "\code{classical}", the estimator is the classical
  instrumental variables estimator based on the sample mean and sample
  variance-covariance matrix.

  If the model contains dummy variables (i.e., \code{dummies != NULL}),
  RIV is computed using an iterative algorithm called "\eqn{L_1}-RIV".
  Briefly, \eqn{L_1}-RIV estimates the coefficients of the dummies using
  an \eqn{L_1}-estimator and the coefficients of the continuous
  covariates using the original RIV. See Cohen Freue et al. for more
  details.
}

\value{
  A list with components:

  \item{Summary.Table}{Matrix of information available about the
    estimator. It contains regression coefficients, and, for
    \code{method} = "\code{S-est}" and "\code{classical}" only, columns
    for the standard error, t-statistic, and p-value.}
  \item{VC}{estimated variance-covariance matrix, computed only if
    \code{method} = "\code{S-est}" or "\code{classical}"}
  \item{MD}{Mahalanobis of each observation to the multivariate location
    S-estimator with respect to the scatter S-estimator (only computed
    if \code{method} = "\code{S-est}").}
  \item{MSE}{vector of three components, computed only if \code{method}
    = "\code{S-est}" or "\code{classical}":
    \itemize{
      \item \code{sigma.hat1}: the mean square error estimation;
      \item \code{sigma.hat2}: the mean square error estimation taking into
        account the weights associated to each observation (only
        computed if \code{method} = "\code{S-est}" and coefficents of
        endogenous variables are exactly identified, i.e., the number of
        instruments is equal to the number of endogenous variables);
      \item \code{sigma.hat3}: the squared median absolute deviation (only
        computed if \code{method} = "\code{S-est}").
    }
  }
  \item{weight}{the weights assigned to each observation (only computed
    if \code{method} = "\code{S-est}" and the number of instruments is
    equal to the number of endogenous variables.}
}

\references{
  LOPUHAA H.P. (1989). \emph{On the Relation between S-estimators and M-estimators of Multivariate Location and Covariance. Ann. Statist.} \bold{17} 1662-1683. 

  COHEN-FREUE, G.V., ORTIZ-MOLINA, H., and ZAMAR, R.H. (2012) \emph{A
    Natural Robustification of the Ordinary Instrumental Variables
    Estimator.} Submitted to \emph{Biometrics}.
}

\seealso{\code{\link[rrcov]{CovSest}}, \code{\link[rrcov]{CovSde}}, \code{\link[rrcov]{CovMcd}}}

\examples{

## load data earthquake: the first column contains the response Y, the
## second the endogenous variable X and the third column is the
## instrument Zinst.
data(earthquake) 
riv.eq <- riv(earthquake[,1],earthquake[,2],NULL,earthquake[,3])

## plot of the RIV estimates and the outlying observations are
## identified by filled points
plot(earthquake[,2],earthquake[,1],xlab="X",ylab="Y",cex=1.5)
abline(riv.eq$Summary.Table[,1])
case <- 1:62
text(earthquake[riv.eq$MD>sqrt(qchisq(0.95, 3)),2],
     earthquake[riv.eq$MD>sqrt(qchisq(0.95, 3)),1],
     case[riv.eq$MD>sqrt(qchisq(0.95, 3))],
     pos=c(1,4,2,2,3,2))
points(earthquake[riv.eq$MD>sqrt(qchisq(0.95, 3)),2],
       earthquake[riv.eq$MD>sqrt(qchisq(0.95, 3)),1],
       cex=1.5,pch=19)

## Mahalanobis distances of the S-estimator of multivariate location and
## covariance (CovSest)
case <- 1:62
plot(riv.eq$MD,xlab="Case number",ylab="Mahalanobis Distance",cex = 1.5)
abline(h=sqrt(qchisq(0.95, 3)))
text(case[riv.eq$MD>sqrt(qchisq(0.95, 3))],
     riv.eq$MD[riv.eq$MD>sqrt(qchisq(0.95, 3))],
     case[riv.eq$MD>sqrt(qchisq(0.95, 3))],
     pos=c(3,3,3,3,2,3))

## Weights associated to each observation 
plot(riv.eq$weight,xlab="Case number",ylab="Weight",cex=1.5)
abline(h=0)
text(case[riv.eq$MD>sqrt(qchisq(0.95, 3))],
     riv.eq$weight[riv.eq$MD>sqrt(qchisq(0.95, 3))],
     case[riv.eq$MD>sqrt(qchisq(0.95, 3))],
     pos=c(3,3,3,2,3,3),
     col=c("black","black","red","black","red","black"))
points(case[riv.eq$weight==0],
       riv.eq$weight[riv.eq$weight==0],
       col="red",cex=1.5,pch=19)


## load data mortality
data(mortality) 
Y <- as.matrix(mortality[,1])            ## M070
Xex <- as.matrix(mortality[,c(2,3,5,6)]) ## MAGE,CI68,DENS,NONW 
Xend <- as.matrix(mortality[,4])         ## MDOC  
colnames(Xend) <- colnames(mortality)[4] 
Zinst <- as.matrix(mortality[,7:8])      ## EDUC,IN69

## Classical instrumental variables estimator
riv(Y, Xend, Xex, Zinst, method="classical")
}

\keyword{robust}
\keyword{multivariate}
