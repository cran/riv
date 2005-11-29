\name{riv}
\alias{riv}
\title{Robust instrumental variables estimator}
\description{
Finds robust instrumental variables estimator using high breakdown point S-estimator of multivariate
location and covariance.
}
\usage{
riv(Y, Xex = FALSE, Xend, W, intercept = TRUE, nsamp = 500, bdp = 0.5, approx = 1000)
}
\arguments{
  \item{Y}{a vector of responses, say of length n.}
  \item{Xex}{a matrix of the exogenous variables, i.e. the variables which are pretented to be uncorrelated with the error. Default = FALSE, i.e. no exogenous
  variable is present in the model.}
  \item{Xend}{a matrix of the endogenous variables, i.e. the variables which are pretended to be correlated with the error.}
  \item{W}{a matrix of instruments. NOTE: we must have an instrument for each endogenous variable.} 
  \item{intercept}{if TRUE (default) a model with intercept is considered.}
  \item{nsamp}{number of subsamples to compute the S-estimator of multivariate location and covariance.}
  \item{bdp}{Breakdown point value of the S-estimator must be 0.15, 0.25 or 0.5 (default).}
  \item{approx}{number of observations generated from the normally distribution with mean zero and the identity matrix as variance
      (default = 1000).}
}
\details{
This estimator is constructed by using a robust S-estimator of multivariate location and covariance based on the Tukey's biweight function (see \code{slc}).

Only the case in which the number of instruments is equal to the number of endogenous variables can be considered with this method, i.e. for each endogenous variable
exact one instrument has to be provided. The matrix W and Xend has to be given appropriate that means the endogenous variable xend\eqn{_i} and the 
associated instrument w\eqn{_i} must be put in the same column respectively.

The estimated variance and covariance matrix is based on the empirical influence function. See references for more details. 
}
\value{
  A list with components
  \item{coefficients}{regression coefficients. if intercept = TRUE the first coefficient is the intercept.}
  \item{VC}{estimated variance-covariance matrix.}
  \item{MD}{Mahalanobis distances from the estimator.}
  \item{MSE}{vector of three components: sigma.hat1 = the mean squared error estimation, 
   sigma.hat2 = the mean squared error estimation taking into account the weights associated to each observation,
   sigma.hat3 = the squared median absolute deviation.}
  \item{weight}{the weights assigned to each observation.}
}
\references{
LOPUHAÄ,H.P. (1989). \emph{On the Relation between S-estimators and M-estimators of Multivariate Location and Covariance. Ann. Statist.} \bold{17} 1662-1683. 

COHEN-FREUE,G.V. and ZAMAR,R.H. (2005). \emph{A Robust Instrumental Variables Estimator.}
}

\author{ 
Kaufmann B. \email{beat.kaufmann@epfl.ch}

Cohen-Freue G.V. \email{gcohen@stat.ubs.ca}
}
\seealso{\code{\link{slc}}}
\examples{
## load data earthquake: the first column contains the response Y, the second the endogenous variable X 
## and the third column is the instrument W.
data(earthquake) 
riv.eq <- riv(earthquake[,1],FALSE,earthquake[,2],earthquake[,3]) 

## plot of the riv estimates and the outlying observations are identified by filled points
plot(earthquake[,2],earthquake[,1],xlab="X",ylab="Y",cex=1.5)
abline(riv.eq$Summary.Table[,1])
case <- 1:62
text(earthquake[riv.eq$MD>sqrt(qchisq(0.95, 3)),2],earthquake[riv.eq$MD>sqrt(qchisq(0.95, 3)),1],
case[riv.eq$MD>sqrt(qchisq(0.95, 3))],pos=c(1,4,2,2,3,2))
points(earthquake[riv.eq$MD>sqrt(qchisq(0.95, 3)),2],earthquake[riv.eq$MD>sqrt(qchisq(0.95, 3)),1],cex=1.5,pch=19)

## Mahalanobis distances of the S-estiamtor of multivariate location and covariance (slc)

case <- 1:62
plot(riv.eq$MD,xlab="Case number",ylab="Mahalanobis Distance",cex = 1.5)
abline(h=sqrt(qchisq(0.95, 3)))
text(case[riv.eq$MD>sqrt(qchisq(0.95, 3))],riv.eq$MD[riv.eq$MD>sqrt(qchisq(0.95, 3))],
case[riv.eq$MD>sqrt(qchisq(0.95, 3))],pos=c(3,3,3,3,2,3))

## Weights associated to each observation 

plot(riv.eq$weight,xlab="Case number",ylab="Weight",cex=1.5)
abline(h=0)
text(case[riv.eq$MD>sqrt(qchisq(0.95, 3))],riv.eq$weight[riv.eq$MD>sqrt(qchisq(0.95, 3))],
case[riv.eq$MD>sqrt(qchisq(0.95, 3))],pos=c(3,3,3,2,3,3),col=c("black","black","red","black","red","black"))
points(case[riv.eq$weight==0],riv.eq$weight[riv.eq$weight==0],col="red",cex=1.5,pch=19)

## load data mortality

data(mortality) 
Y <- as.matrix(mortality[,1])            ## M070
Xex <- as.matrix(mortality[,c(2,3,5,6)]) ## MAGE,CI68,DENS,NONW 
Xend <- as.matrix(mortality[,4])         ## MDOC  
colnames(Xend) <- colnames(mortality)[4] 
W <- as.matrix(mortality[,8])            ## IN69
riv(Y,Xex,Xend,W)
}

\keyword{robust}
\keyword{multivariate}
