"sestck" <-
function(x,start,c,k,tol)
{
if (start>0) s<-start
else  
  {
  x <- as.matrix(x)
  x <- abs(x)
  n<-nrow(x)
  p<-ncol(x)
  x<-apply(x,2,sort)       
  if (floor(n/2)==n/2) s<-(x[n/2,]+x[(n+2)/2,])/2
  else s<-x[(n+1)/2,]
  s <- s/0.6745
  }
crit<-2*tol
rhoold<-mean(rho.biweight(x/s,c))-k
while (crit>=tol){
      delta<-rhoold/mean(psi.bisquare(x/s,c)*(x^2/(s^3)))
      isqu<-1
      okay<-0
      while((isqu<10)&(okay!=1)){
          rhonew<-mean(rho.biweight(x/(s+delta),c))-k
          if (abs(rhonew)<abs(rhoold)){
              s<-s+delta
              okay<-1
                           } 
          else {
              delta<-delta/2 
              isqu<-isqu+1
               }
       }
     if (isqu==10) crit<-0
     else crit<-(abs(rhoold)-abs(rhonew))/max(abs(rhonew),tol)
     rhoold<-rhonew
     }
scale<-abs(s)
scale
}

