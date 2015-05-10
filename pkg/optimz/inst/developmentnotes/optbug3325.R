## optbug3325.R
rm(list=ls())
require(optimz)

sn.dev <- function(cp, X, y, trace=FALSE)
{ # -2*logL for centred parameters  
  m <- ncol(X)
  if(abs(cp[length(cp)])> 0.99527) {print(cp); stop("cp outside bounds")}
  dp <- as.vector(cp.to.dp(cp))
  location <- X %*% as.matrix(dp[1:m])
  scale <- dp[m+1]
  # AVOID: logL <- sum(log(dsn(y,location,dp[m+1],dp[m+2])))
  z <- (y-location)/scale
  nlogL <- (length(y)*log(2.506628274631*scale) + 0.5*sum(z^2)
            - sum(zeta(0,dp[m+2]*z)))
  if(trace) {cat("sn.dev: (cp,dev)="); print(c(cp,2*nlogL))}
  return(2*nlogL) 
}

sn.dev.gh <- function(cp, X, y, trace=FALSE, hessian=FALSE)
{
  # computes gradient and hessian of dev=-2*logL for centred parameters 
  # (and observed information matrix);
  m  <- ncol(X)
  n  <- nrow(X)
  np <- m+2
  score <- rep(NA,np)
  info  <- matrix(NA,np,np)
  beta <- cp[1:m]
  sigma <- cp[m+1]
  gamma1 <- cp[m+2]
  lambda <- gamma1.to.lambda(gamma1)
  # dp<-cp.to.dp(c(beta,sigma,gamma1))
  # info.dp <- sn.info(dp,y)$info.dp
  mu <- as.vector(X %*% as.matrix(beta))
  d  <- y-mu
  r  <- d/sigma
  E.Z<- lambda*sqrt(2/(pi*(1+lambda^2)))
  s.Z<- sqrt(1-E.Z^2)
  z  <- E.Z+s.Z*r
  p1 <- as.vector(zeta(1,lambda*z))
  p2 <- as.vector(zeta(2,lambda*z))
  omega<- sigma/s.Z
  w    <- lambda*p1-E.Z
  DE.Z <- sqrt(2/pi)/(1+lambda^2)^1.5
  Ds.Z <- (-E.Z/s.Z)*DE.Z
  Dz   <- DE.Z + r*Ds.Z
  DDE.Z<- (-3)*E.Z/(1+lambda^2)^2
  DDs.Z<- -((DE.Z*s.Z-E.Z*Ds.Z)*DE.Z/s.Z^2+E.Z*DDE.Z/s.Z)
  DDz  <- DDE.Z + r*DDs.Z
  score[1:m] <- omega^(-2)*t(X) %*% as.matrix(y-mu-omega*w) 
  score[m+1] <- (-n)/sigma+s.Z*sum(d*(z-p1*lambda))/sigma^2
  score[m+2] <- score.l <- n*Ds.Z/s.Z-sum(z*Dz)+sum(p1*(z+lambda*Dz))
  Dg.Dl <-1.5*(4-pi)*E.Z^2*(DE.Z*s.Z-E.Z*Ds.Z)/s.Z^4
  R <- E.Z/s.Z
  T <- sqrt(2/pi-(1-2/pi)*R^2)
  Dl.Dg <- 2*(T/(T*R)^2+(1-2/pi)/T^3)/(3*(4-pi))
  R. <- 2/(3*R^2 * (4-pi))
  T. <- (-R)*R.*(1-2/pi)/T
  DDl.Dg <- (-2/(3*(4-pi))) * (T./(R*T)^2+2*R./(T*R^3)+3*(1-2/pi)*T./T^4)
  score[m+2] <- score[m+2]/Dg.Dl  # convert deriv wrt lamda to gamma1 
  gradient <- (-2)*score
  if(hessian){
     info[1:m,1:m] <- omega^(-2) * t(X) %*% ((1-lambda^2*p2)*X)
     info[1:m,m+1] <- info[m+1,1:m] <- 
            s.Z* t(X) %*% as.matrix((z-lambda*p1)+d*(1-lambda^2*p2)*
            s.Z/sigma)/sigma^2
     info[m+1,m+1] <- (-n)/sigma^2+2*s.Z*sum(d*(z-lambda*p1))/sigma^3 +
            s.Z^2*sum(d*(1-lambda^2*p2)*d)/sigma^4
     info[1:m,m+2] <- info[m+2,1:m] <- 
            t(X)%*%(-2*Ds.Z*d/omega+Ds.Z*w+s.Z*(p1+lambda*p2*(z+lambda*Dz)
            -DE.Z))/sigma 
     info[m+1,m+2] <- info[m+2,m+1] <- 
            -sum(d*(Ds.Z*(z-lambda*p1)+s.Z*(Dz-p1-p2*lambda*(z+lambda*Dz))
             ))/sigma^2
     info[m+2,m+2] <- (n*(-DDs.Z*s.Z+Ds.Z^2)/s.Z^2+sum(Dz^2+z*DDz)-
            sum(p2*(z+lambda*Dz)^2)- sum(p1*(2*Dz+lambda*DDz)))
     info[np,] <- info[np,]/Dg.Dl # convert info wrt lamda to gamma1 
     info[,np] <- info[,np]*Dl.Dg # an equivalent form of the above
     info[np,np] <- info[np,np]-score.l*DDl.Dg
     }
  attr(gradient,"hessian") <- 2*info
  if(trace) {cat("sn.dev.gh: gradient="); print(-2*score)}
  return(gradient)
}

cp.to.dp <- function(param){
  # converts centred parameters cp=(mu,sigma,gamma1)
  # to direct parameters dp=(xi,omega,lambda)
  # Note:  mu can be m-dimensional, the other must be scalars
  b <- sqrt(2/pi)
  m <- length(param)-2
  gamma1 <- param[m+2]
  if(abs(gamma1)>0.9952719) stop("abs(gamma1)>0.9952719 ")
  A <- sign(gamma1)*(abs(2*gamma1/(4-pi)))^(1/3)
  delta <- A/(b*sqrt(1+A^2))
  lambda <- delta/sqrt(1-delta^2)
  E.Z  <- b*delta
  sd.Z <- sqrt(1-E.Z^2)
  location    <- param[1:m]
  location[1] <- param[1]-param[m+1]*E.Z/sd.Z
  scale <- param[m+1]/sd.Z
  dp    <- c(location,scale,lambda)
  names(dp)[(m+1):(m+2)] <- c("scale","shape")
  if(m==1)  names(dp)[1] <- "location"
  dp
  }

zeta <- function(k,x){# k integer \in (0,4)
  k <- as.integer(k)
  na <- is.na(x)
  x <- replace(x,na,0)
  if(any(abs(x)==Inf)) stop("Inf not allowed")
  # funzionerebbe per k=0 e 1, ma non per k>1
  ok <- (-35<x)
  if(k>2 & sum(!ok)>0) na<- (na | x < (-35))
  if(k==0)  
    {ax <- (-x[!ok])
    ay  <- (-0.918938533204673)-0.5*ax^2-log(ax)
    y   <- rep(NA,length(x))
    y[ok] <- log(2*pnorm(x[ok]))
    y[!ok]<- ay
    }
  else {if(k==1) {y  <- (-x)*(1+1/x^2)
          y[ok]<-dnorm(x[ok])/pnorm(x[ok]) }
    else { if(k==2)  y<-(-zeta(1,x)*(x+zeta(1,x)))
      else{ if(k==3)  y<-(-zeta(2,x)*(x+zeta(1,x))-zeta(1,x)*(1+zeta(2,x)))
        else{ if(k==4)  
           y<-(-zeta(3,x)*(x+2*zeta(1,x))-2*zeta(2,x)*(1+zeta(2,x)))
        else stop("k must be integer in (0,4)") }}}}
  replace(y,na,NA)
}

gamma1.to.lambda<- function(gamma1){
  max.gamma1 <- 0.5*(4-pi)*(2/(pi-2))^1.5
  na <- (abs(gamma1)>max.gamma1)
  if(any(na)) warning("NAs generated") 
  gamma1<-replace(gamma1,na,NA)
  a    <- sign(gamma1)*(2*abs(gamma1)/(4-pi))^0.33333
  delta<- sqrt(pi/2)*a/sqrt(1+a^2)
  lambda<-delta/sqrt(1-delta^2)
  as.vector(lambda)
}

monica2 <-  structure(c(-1.13886379255452, -1.49413719201167, -1.87841064955904, 
-0.513738594928947, -1.00896950224601, -1.6532307897784, -1.30938967543405, 
-1.06041025293133, -1.81471688461199, -1.33766295196981, -1.64330728995432, 
-1.06469899752166, -1.00001063084463, -1.39630202884595, -0.69585045438959, 
-1.69203796703315, -1.52865988854611, -1.62709073708234, -0.934436685805195, 
-1.11123590657362, -1.67332502980791, -1.16134216754217, -1.27483304830845, 
-1.83766242751329, -1.56429386551254, -0.70663591112117, -0.806879736670856, 
-1.39706539761330, -1.48444503984765, -1.48192041745309, -0.763350185813698, 
-1.42340943250413, -1.5034850496403, -0.916057843936642, -0.797650717731944, 
-1.90218776624194, -1.89573584905280, -1.55031254908193, -1.55488970245028, 
-1.17137569123985, -1.81490209565657, -1.82974300222698, -0.797407302025425, 
-1.16514651844490, -1.31007112485982, -1.71325589588157, -1.54253707834837, 
-1.28190616471927, -1.28990272862791, -0.593534374170964),
          parameters = c(-1.71631422551004, 0.5, 2))

#----
# cp <- c(-1.322168500,  0.373810085,  0.326283975)
y <- monica2
X <-  matrix(1, nrow=length(monica2))
n <- length(y)
m <- ncol(X)
    qrX <- qr(X)
    s <- sqrt(sum(qr.resid(qrX, y)^2)/n)
    gamma1 <- sum(qr.resid(qrX, y)^3)/(n*s^3)
    if(abs(gamma1)>0.99527) gamma1<- sign(gamma1)*0.95
    cp <- c(qr.coef(qrX,y), s, gamma1)
opt<- optimx(cp, fn=sn.dev, gr=sn.dev.gh, method=c("L-BFGS-B","lbfgsb3"),
         lower=c(-rep(Inf,m),1e-10, -0.99527), 
         upper=c(rep(Inf,m), Inf, 0.99527), 
         control=list(iter.max=100, abs.tol=1e-5, trace=1),
         X=X, y=y, hessian=FALSE)
opt
tmp <- readline("Now try all meths")

optall<- optimx(cp, fn=sn.dev, gr=sn.dev.gh, method="all",
         lower=c(-rep(Inf,m),1e-10, -0.99527), 
         upper=c(rep(Inf,m), Inf, 0.99527), 
         control=list(iter.max=100, abs.tol=1e-5, trace=1),
         X=X, y=y, hessian=FALSE)
summary(optall, order=value)
