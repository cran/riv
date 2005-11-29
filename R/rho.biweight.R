"rho.biweight" <-
function(x,c){
   hulp<-x^2/2-x^4/(2*c^2)+x^6/(6*c^4)
   rho<-hulp*(abs(x)<c)+c^2/6*(abs(x)>=c)
   rho}

