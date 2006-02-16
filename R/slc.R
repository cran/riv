"slc" <-
function(x,nsamp=500,bdp=0.5)
{
Tbsb<-function(c,p){
ksiint<-function(cc,ss,pp){(2^ss)*gamma(ss+pp/2)*pgamma(cc^2/2,ss+pp/2)/gamma(pp/2)}
y1<-ksiint(c,1,p)*3/c-ksiint(c,2,p)*3/(c^3)+ksiint(c,3,p)/(c^5)
y2<-c*(1-pchisq(c^2,p))
res<-y1+y2
res 
}

if (!require(MASS))
stop("cannot load required library MASS")

if (is.null(dim(x))) x=as.matrix(x)                                                         
if (any(is.na(x))) stop("missing values are not allowed")                                   
if (bdp!=0.5 & bdp != 0.15 & bdp != 0.25) stop("bdp must be 0.15,0.25 or 0.5")              
n<-nrow(x)
p<-ncol(x)
if (p > n) stop("number of observations must be greater than number of variables")          
tol<-10^(-5)
s<-10^(11)

tbdp<- sqrt(qchisq(1-bdp,p))
maxit<-1000
eps<-10^(-8)
diff<-10^6
ctest<-tbdp
iter<-1
while ((diff>eps)& (iter<maxit)) 
  {
  cold<-ctest
  ctest<-Tbsb(cold,p)/bdp
  diff<-abs(cold-ctest)
  eter<-iter+1
  }
c<-ctest

k<-(c/6)*Tbsb(c,p)
la<-1

for (i in 1:nsamp)
  {    # global improvement
  ranset<- sample(1:n,p+1)                                                                  
  xj<-as.matrix(x[ranset,])                                                                 
  mu<-apply(xj,2,mean)
  cov<-var(xj)*(nrow(xj)-1)/nrow(xj)
  determ<-det(cov)
  if ((determ>10^(-15))&(determ^(1/p)>10^(-5))) 
    {
    cov<-determ^(-1/p)*cov
    if (i>ceiling(nsamp/5)) 
      {
      if (i==ceiling(nsamp/2)) la<-2
      if (i==ceiling(nsamp*.8)) la<-4
      random<- runif(1) 
      random<-random^la
      mu<-random*mu+(1-random)*muopt                                                      
      cov<-random*cov+(1-random)*covopt
      determ<-det(cov)
      cov<-determ^(-1/p)*cov
      }
    md<-mahalanobis(x,mu,cov,inverted = FALSE, tol.inv =.Machine$double.eps)
    md<-md^(1/2)
    if (mean(rho.biweight(md/s,c))<k) 
      {
      if (s<5*10^10) s<-sestck(md,s,c,k,tol)
      else s<-sestck(md,0,c,k,tol)
      muopt<-mu
      covopt<-cov
      mdopt<-md                                                                            
      psi<-psi.bisquare(md,s*c)*md
      u<-psi.bisquare(md,s*c)
      ubig<-matrix(t(u),nrow=length(u),ncol=p,byrow=FALSE)
      aux<-(ubig*x)/mean(u)
      mu<-apply(aux,2,mean)
      xcenter<-t(t(x)-mu)
      cov<-t(ubig*xcenter)%*%xcenter
      cov<-det(cov)^(-1/p)*cov
      okay<-0
      jj<-1
      while ((jj<3)&(okay!=1)) 
        {
        jj<-jj+1
        md<-mahalanobis(x,mu,cov,tol.inv =.Machine$double.eps)
        md<-md^(1/2)
        if (mean(rho.biweight(md/s,c))<k)
          {
          muopt<-mu
          covopt <-cov
          mdopt <-md
          okay<-1
          if (s<5*10^10) s<-sestck(md,s,c,k,tol)
          else s<-sestck(md,0,c,k,tol)
          }
        else 
          {
          mu<-(mu+muopt)/2
          cov<-determ^(-1/p)*(cov+covopt)/2
          }
        }
      }
    } 
  }
res.mean<-muopt
res.covariance<-covopt*s^2
res.distances<-mdopt/s
res.scale<-s
list(location=res.mean,covariance=res.covariance,distances=res.distances,scale=res.scale,c=c)
}

