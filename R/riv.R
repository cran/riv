"riv" <-
function(Y,Xex=FALSE,Xend,W,intercept=TRUE,method=c("robust","classical"),nsamp=500,bdp=.5,approx=1000)
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

  multarray <- function(x,B)
  {
  result <- array(NA,dim=c(dim(x)[1],ncol(B),dim(x)[3]))
  for (i in 1 : dim(x)[3]) result[,,i] <- x[,,i] %*% B
  result
  }
 
  multarray2 <- function(x,B)
  {
  result <- array(NA,dim=c(nrow(B),dim(x)[2],dim(x)[3]))
  for (i in 1 : dim(x)[3]) result[,,i] <- B %*% x[,,i]
  result
  }

  if (!require(MASS))
  stop("cannot load required library MASS")
  n=length(Y)
  Xend=as.matrix(Xend)
  W=as.matrix(W)
  kend = ncol(Xend)
  k=ncol(W)
  if (kend > k) 
  stop("the number of instruments must be equal or bigger to the number of endogenous variables")
  if (is.logical(Xex))  
  {      
    X = Xend
    if(is.null(colnames(X))) colnames(X) <- paste("X",seq(ncol(X)),sep="")
    Z = cbind(Xend,W,Y)
  }else
  {
    Xex=as.matrix(Xex)
    X = cbind(Xend,Xex) 
    if(is.null(colnames(X))) 
    colnames(X) <- c(paste("Xend",seq(ncol(Xend)),sep=""),paste("Xex",seq(ncol(Xex)),sep=""))
    else
    {
      if(is.null(colnames(Xend))) colnames(X) <- c(paste("Xend",seq(ncol(Xend)),sep=""),colnames(Xex))
      if(is.null(colnames(Xex))) colnames(X) <- c(colnames(Xend),paste("Xex",seq(ncol(Xex)),sep=""))
    }
    Z = cbind(Xend,W,Xex,Y)
  }
  p=ncol(X)
  r=k+p+1
  if (any(is.na(Z))) stop("missing values are not allowed")   
  if (any(is.na(X))) stop("missing values are not allowed")                               
  if (bdp!=0.5 & bdp != 0.15 & bdp != 0.25) stop("bdp must be 0.15,0.25 or 0.5")             
  method <- match.arg(method)
  
  if (method == "robust")
  {
    res1=slc(Z,nsamp,bdp)                                                                       
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
 
    if (k == kend)
    {
      hulp=.5*MD^2-(1.5*MD^4)/c^2+(5*MD^6)/(6*c^4)+Econst$b0
      v=(Econst$b0-(c^2)/6)+(hulp-(Econst$b0-(c^2)/6))*(abs(MD)<c)
      weight=r*psi.bisquare(MD,c)/mean(v)
      weight=weight/sum(weight)
    }    

    B<- t(chol(V))
    Binv=solve(B) 
    r=ncol(Z)
    z0=Binv%*%(t(Z)-L)                                                                       
    norm.z0=sqrt(apply(z0^2,2,sum))                                                          
    IF.M0=t(z0)*(1/Econst$nu)*psi.bisquare(norm.z0,c)
    IFL=B%*%t(IF.M0)
    IFV=as.vector(apply(z0,2,IF.obs,B,Econst,c))
    dim(IFV)=c(nrow(z0),nrow(z0),ncol(z0))
    
    Vm=matrix(V[(kend+1):(nrow(V)-1),-((kend+1):(kend+k))],nrow=nrow(V)-kend-1)     
    Swx=Vm[,-ncol(Vm)]
    Sxw=t(Swx)
    Sww=matrix(V[(kend+1):(nrow(V)-1),(kend+1):(nrow(V)-1)],nrow=nrow(V)-kend-1)
    Swy=Vm[,ncol(Vm)]
    Lm=L[-((kend+1):(kend+k))]
    Mx=Lm[1:(length(Lm)-1)]
    My=Lm[length(Lm)]
    part1 <- Sxw%*%solve(Sww)%*%Swx
    b1<- solve(part1) %*% Sxw %*% solve(Sww) %*% Swy
    biv=b1
    resid=Y-X%*%biv

    if (intercept)                             
    {
      b0<-My-sum(b1*Mx)
      biv<-matrix(rbind(b0,b1),ncol=1,dimnames=NULL)
      resid=Y-b0-X%*%b1
    }
    df.riv= n-p-intercept
    IFL=IFL[-((k+1):(kend+k)),]
    IF.Mx=IFL[1:(length(Lm)-1),]
    IF.My=IFL[length(Lm),]

    Swwinv <- solve(Sww)
    A = Sxw %*% Swwinv %*% Swx
    B = Sxw %*% Swwinv %*% Swy

    if (p+k>2) 
    {
      IF.Swx=IFV[(kend+1):(r-1),-c((kend+1):(kend+k),r),]
      IF.Sxw=IFV[-c((kend+1):(kend+k),r),(kend+1):(r-1),]
      IF.Sww=IFV[(kend+1):(r-1),(kend+1):(r-1),]
      IF.Swy=IFV[(kend+1):(r-1),r,]
      IF.prov <- array(NA,dim=c(p-kend+k,1,n))    
      for (i in 1:n) IF.prov[,,i] <- t(IF.Swy[,i])
      IF.Swy <- IF.prov

      if (p==1)
      { 
      IF.Xprov <- array(NA,dim=c(p,k,n))
      IF.Wprov <- array(NA,dim=c(k,p,n))
      for (i in 1:n) 
        { 
        IF.Xprov[,,i] <- t(IF.Sxw[,i])
        IF.Wprov[,,i] <- IF.Swx[,i]
        }   
      IF.Swx <- IF.Wprov
      IF.Sxw <- IF.Xprov
      }

      IF.A1 <- multarray(IF.Sxw,solve(Sww)%*%Swx)
      IF.A2a <- multarray(IF.Sww,solve(Sww)%*%Swx)
      IF.A2 <- -multarray2(IF.A2a,Sxw%*%solve(Sww))
    
      IF.A3 <- IF.A1
      for (i in 1:n) IF.A3[,,i] <- t(IF.A1[,,i])
   
      IF.A <- IF.A1 + IF.A2  + IF.A3
      IF.solveA <- multarray(IF.A,solve(A)%*%B)
      IF.solveA <- -multarray2(IF.solveA,solve(A))
    
      IF.B1 <- multarray(IF.Sxw,solve(Sww)%*%Swy)
      IF.B2a <- multarray(IF.Sww,solve(Sww)%*%Swy)
      IF.B2 <- -multarray2(IF.B2a,Sxw%*%solve(Sww))
      IF.B3 <- multarray2(IF.Swy,Sxw%*%solve(Sww))
     
      IF.B <- IF.B1 + IF.B2 + IF.B3
      IF.B <- multarray2(IF.B,solve(A))

      IF.b1.array <- IF.solveA + IF.B
   
      IF.b1 <- t(matrix(IF.b1.array,ncol=p,byrow=TRUE))
      IFriv= IF.b1
      if (intercept)
      {
        term1 <- multarray2(IF.b1.array,t(Mx))
        term1 <- c(term1[1,1,])
        term2 <- as.vector(t(b1)%*%IF.Mx)
        IF.b0 <- IF.My - term1 - term2
        IF.b0 <- matrix(IF.b0,nrow=1)
        IFriv <- rbind(IF.b0,IF.b1)
      }
    }  
    else  
    {
      IF.Swx=IFV[(kend+1):(r-1),-c((kend+1):(kend+k),r),]
      IF.Sxw=IFV[-c((kend+1):(kend+k),r),(kend+1):(r-1),]
      IF.Sww=IFV[(kend+1):(kend+k),(kend+1):(kend+k),]
      IF.Swy=IFV[(kend+1):(kend+k),r,]
   
      IF.A1 <- IF.Sxw*solve(Sww)*Swx
      IF.A2 <- -Sxw*solve(Sww)*IF.Sww*solve(Sww)*Swx
      IF.A3 <- Sxw*solve(Sww)*IF.Swx
      IF.A <- IF.A1 + IF.A2  + IF.A3
      IF.A <- -solve(A)*IF.A*solve(A)*B
  
      IF.B1 <- IF.Sxw*solve(Sww)*Swy 
      IF.B2 <- -Sxw*solve(Sww)*IF.Sww*solve(Sww)*Swy 
      IF.B3 <- Sxw*solve(Sww)*IF.Swy

      IF.B <- IF.B1 + IF.B2 + IF.B3
      IF.B <- solve(A)*IF.B
    
      IF.b1 <- IF.A + IF.B
      Ifriv <- IF.b1
      if (intercept)
      {   
        IF.b0 <- IF.My - IF.b1*Mx - b1*IF.Mx
        IFriv <- rbind(IF.b0,IF.b1)
      }
    }
    n=ncol(IFriv)
    A=as.vector(apply(IFriv,2,function(x) {x%*%t(x)}))
    dim(A)=c((p+intercept),(p+intercept),n)
    res5=apply(A,c(1,2),mean)/n
    if (intercept) dimnames(res5) = list(c("intercept",colnames(X)),c("intercept",colnames(X)))
    else dimnames(res5) = list(colnames(X),colnames(X))
    sd.riv=sqrt(diag(res5))
    tval=biv/sd.riv
    pv=2*(1-pt(abs(tval),df.riv))
    sigma.hat1=as.numeric(crossprod(resid)/df.riv)
    if (k == kend)
    {
    w.resid=resid*weight
    sigma.hat2=as.numeric(crossprod(w.resid)/df.riv)
    }
    sigma.hat3=(mad(resid))^2
    tabRIV=cbind(biv,sd.riv,tval,pv)
    colnames(tabRIV)=c("Coef","Std.Err.","t","p.values")
    if (k == kend) res=list(Summary.Table=tabRIV,VC=res5,MSE=c(sigma.hat1,sigma.hat2,sigma.hat3),MD=MD,weight=weight)
    else res=list(Summary.Table=tabRIV,VC=res5,MSE=c(sigma.hat1,sigma.hat3),MD=MD)
    res
  }
  else if (method == "classical")
  {
    if (is.logical(Xex)) 
    { 			
    X = Xend
    if(is.null(colnames(X))) colnames(X) <- paste("Xend",seq(ncol(X)),sep="")
    Xnam=X
    Z = W
  }
  else  
  {
    Xex=as.matrix(Xex)
    X = cbind(Xend,Xex)
    if(is.null(colnames(X))) 
    colnames(X) <- c(paste("Xend",seq(ncol(Xend)),sep=""),paste("Xex",seq(ncol(Xex)),sep=""))
    else
    {
      if(is.null(colnames(Xend))) colnames(X) <- c(paste("Xend",seq(ncol(Xend)),sep=""),colnames(Xex))
      if(is.null(colnames(Xex))) colnames(X) <- c(colnames(Xend),paste("Xex",seq(ncol(Xex)),sep=""))
    }
    Z = cbind(W,Xex)
  }
  p=ncol(X)
    df.oiv=n-p-intercept
    if (intercept)
    { 
      Z=cbind(1,Z)
      x=cbind(1,X)
      zinv=solve(crossprod(Z))
      Pz=Z%*%zinv%*%t(Z)
      XPz=t(x)%*%Pz%*%x
      inv.XPz=solve(XPz)
      XPy=t(x)%*%Pz%*%Y
      beta=inv.XPz%*%XPy
      resid=Y-x%*%beta
    }
    else
    {
      zinv=solve(crossprod(Z))
      Pz=Z%*%zinv%*%t(Z)
      XPz=t(X)%*%Pz%*%X
      inv.XPz=solve(XPz)
      XPy=t(X)%*%Pz%*%Y
      beta=inv.XPz%*%XPy
      resid=Y-X%*%beta
    }
    sigma.hat=as.numeric(crossprod(resid)/df.oiv)
    var.oiv=sigma.hat*inv.XPz
    beta.oiv=beta
    sd.oiv=sqrt(diag(var.oiv))
    tval=beta.oiv/sd.oiv
    pv=2*(1-pt(abs(tval),df.oiv))
  
    tabOIV=cbind(beta.oiv,sd.oiv,tval,pv)
    colnames(tabOIV)=c("Coef","Std.Err.","t","p.values")
    if (intercept)
    rownames(tabOIV) <- colnames(var.oiv) <- rownames(var.oiv) <- c("Intercept",colnames(X))
    else rownames(tabOIV) <- colnames(var.oiv) <- rownames(var.oiv) <- colnames(X)
    res=return(list(Summary.Table=tabOIV,VC=var.oiv,MSE=sigma.hat))
    res
  }
  else stop("'method' is unknown") 
}
