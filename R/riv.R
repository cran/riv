"riv" <-
function(Y,Xex=FALSE,Xend,W,intercept=TRUE,nsamp=500,bdp=.5,approx=1000)
{
  psi.prime=function(x,c)
  {
  hulp=1-6*x^2/(c^2)+5*x^4/(c^4)
  psipr<-hulp*(abs(x)<c)
  psipr
  }

  IF.obs=function(x,B,const,c)
  {
    p=length(x)
    norm.x=sqrt(sum(x^2))
    ifs0= ((x%*%t(x)/norm.x^2)-diag(p)/p)*psi.bisquare(norm.x,c)*norm.x^2*p/const$gamma1+
    (2/const$gamma3)*(rho.biweight(norm.x,c)-const$b0)*diag(p)
    ifs=B%*%ifs0%*%t(B)
    ifs
  }

  postmult=function(x,B){x%*%B}    

  mult = function(x) {x%*%t(x)}

  if (!require(MASS))
  stop("cannot load required library MASS")
  n=length(Y)
  Xend=as.matrix(Xend)
  W=as.matrix(W)
  kend = ncol(Xend)
  k=ncol(W)
  if (kend != k) stop("the number of instruments must be equal to the number of endogenous variables")
  if (is.logical(Xex)) 
  {      
    X = Xend
    if(is.null(colnames(X))) colnames(X) <- paste("X",seq(ncol(X)),sep="")
    z = cbind(Xend,W,Y)
  }
  else  
  {
    Xex=as.matrix(Xex)
    X = cbind(Xex,Xend)
    if(is.null(colnames(Xend))) colnames(X) <- c(paste("Xend",seq(ncol(Xend)),sep=""),colnames(Xex))
    if(is.null(colnames(Xex))) colnames(X) <- c(colnames(Xend),paste("Xex",seq(ncol(Xex)),sep=""))
    if(is.null(colnames(X))) colnames(X) <- paste("X",seq(ncol(X)),sep="")
    z = cbind(Xend,W,Xex,Y)
  }
  p=ncol(X)
  r=k+p+1
  if (any(is.na(z))) stop("missing values are not allowed")                                  
  if (bdp!=0.5 & bdp != 0.15 & bdp != 0.25) stop("bdp must be 0.15,0.25 or 0.5")             

  res1=slc(z,nsamp,bdp)                                                                       
  L=res1[[1]]                                                                                 
  V=res1[[2]]                                                                                 
  MD=res1[[3]]                                                                                
  c = res1[[5]]                                                                               

  x.sample= mvrnorm(approx,rep(0,r),diag(r))
  x.sample= x.sample^2
  x.norms=sqrt(apply(x.sample,1,sum))
  fnu=(1-1/r)*psi.bisquare(x.norms,c)+(1/r)*psi.prime(x.norms,c)
  nu=mean(fnu)
  b=mean(rho.biweight(x.norms,c))
  fg3=psi.bisquare(x.norms,c)*x.norms^2
  g3=mean(fg3)
  fg1=psi.prime(x.norms,c)*x.norms^2+(r+1)*psi.bisquare(x.norms,c)*x.norms^2
  g1=mean(fg1)/(r+2)
  fg2=2*psi.prime(x.norms,c)*x.norms^2+r*psi.bisquare(x.norms,c)*x.norms^2
  g2=mean(fg2)/(2*r*(r+2))
  Econst=list(nu=nu,b0=b,gamma1=g1,gamma2=g2,gamma3=g3)
 
  hulp=.5*MD^2-(1.5*MD^4)/c^2+(5*MD^6)/(6*c^4)+Econst$b0
  v=(Econst$b0-(c^2)/6)+(hulp-(Econst$b0-(c^2)/6))*(abs(MD)<c)
  weight=r*psi.bisquare(MD,c)/mean(v)
  weight=weight/sum(weight)
    
  B<- t(chol(V))                                                                            
  Binv=solve(B)          
  r=ncol(z)
  z0=Binv%*%(t(z)-L)                                                                       
  norm.z0=sqrt(apply(z0^2,2,sum))                                                          
  IF.M0=t(z0)*(1/Econst$nu)*psi.bisquare(norm.z0,c)
  IFL=B%*%t(IF.M0)
  IFV=as.vector(apply(z0,2,IF.obs,B,Econst,c))
  dim(IFV)=c(nrow(z0),nrow(z0),ncol(z0))
    
  Vm=matrix(V[(k+1):(nrow(V)-1),-((k+1):(2*k))],nrow=nrow(V)-k-1)     
  Swx=Vm[,-ncol(Vm)]
  Swxinv=solve(Swx)
  Swy=Vm[,ncol(Vm)]
  Lm=L[-((k+1):(2*k))]
  Mx=Lm[1:(length(Lm)-1)]
  My=Lm[length(Lm)]

  b1<-Swxinv%*%Swy                                                                            
  biv=b1
  df.riv=n-p
  resid=Y-X%*%biv

  if (intercept)                             
  {
    b0<-My-sum(b1*Mx)
    biv<-matrix(rbind(b0,b1),ncol=1,dimnames=NULL)
    resid=Y-b0-X%*%b1
    df.riv=n-p-1
  }
  
  IFL=IFL[-((k+1):(2*k)),]
  IF.Mx=IFL[1:(length(Lm)-1),]
  IF.My=IFL[length(Lm),]

  if (p>1) 
  {
    IFV=IFV[(k+1):(r-1),-((k+1):(2*k)),]
    IF.Swx=IFV[,-ncol(Vm),]
    IF.Swy=IFV[,ncol(Vm),]
    IF.b1=-Swxinv%*%apply(IF.Swx,c(1,3),postmult,b1)+Swxinv%*%IF.Swy
    IF.b0=IF.My-as.vector(t(IF.b1)%*%Mx)-as.vector(t(b1)%*%IF.Mx)
    IFriv=rbind(t(IF.b0),IF.b1)
  }
  else 
  {
    IF.Swx=IFV[(p+1):(2*p),1:p,]
    IF.Swy=IFV[(p+1):(2*p),r,]
    IF.b1=(-IF.Swx*b1+IF.Swy)/Swx
    IF.b0=IF.My-IF.b1*Mx-b1*IF.Mx
    IFriv=rbind(IF.b0,IF.b1) 
  }

  p=nrow(IFriv)-1
  n=ncol(IFriv)
  A=as.vector(apply(IFriv,2,mult))
  dim(A)=c((p+1),(p+1),n)
  res5=apply(A,c(1,2),mean)/n
  if (intercept) dimnames(res5) = list(c("intercept",colnames(X)),c("intercept",colnames(X)))
  else dimnames(res5) = list(colnames(X),colnames(X))
  
  sd.riv=sqrt(diag(res5))
  tval=biv/sd.riv
  pv=2*(1-pt(abs(tval),df.riv))
  sigma.hat1=as.numeric(crossprod(resid)/df.riv)
  w.resid=resid*weight
  sigma.hat2=as.numeric(crossprod(w.resid)/df.riv)
  sigma.hat3=(mad(resid))^2

  tabRIV=cbind(biv,sd.riv,tval,pv)
  colnames(tabRIV)=c("Coef","Std.Err.","t","p.values")

  res=return(list(Summary.Table=tabRIV,VC=res5,MSE=c(sigma.hat1,sigma.hat2,
  sigma.hat3),MD=MD,weight=weight))
  res
}

