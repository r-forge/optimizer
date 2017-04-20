## Optimization test function HOBBS
## ?? refs (put in .doc??)
## Nash and Walker-Smith (1987, 1989) ...


hobbs.f<- function(x){ # # Hobbs weeds problem -- function
    if (abs(12*x[3]) > 500) { # check computability
       fbad<-.Machine$double.xmax
       return(fbad)
    }
    res<-hobbs.res(x)
    f<-sum(res*res)
}


hobbs.res<-function(x){ # Hobbs weeds problem -- residual
# This variant uses looping
    if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
    y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 50.156, 62.948,
         75.995, 91.972)
    t<-1:12
    if(abs(12*x[3])>50) {
       res<-rep(Inf,12)
    } else {
       res<-x[1]/(1+x[2]*exp(-x[3]*t)) - y
    }
}

hobbs.jac<-function(x){ # Jacobian of Hobbs weeds problem
   jj<-matrix(0.0, 12, 3)
   t<-1:12
    yy<-exp(-x[3]*t)
    zz<-1.0/(1+x[2]*yy)
     jj[t,1] <- zz
     jj[t,2] <- -x[1]*zz*zz*yy
     jj[t,3] <- x[1]*zz*zz*yy*x[2]*t
   return(jj)
}

hobbs.g<-function(x){ # gradient of Hobbs weeds problem
    # NOT EFFICIENT TO CALL AGAIN
    jj<-hobbs.jac(x)
    res<-hobbs.res(x)
    gg<-as.vector(2.*t(jj) %*% res)
    return(gg)
}


hobbs.rsd<-function(x) { # Jacobian second derivative
    rsd<-array(0.0, c(12,3,3))
    t<-1:12
    yy<-exp(-x[3]*t)
    zz<-1.0/(1+x[2]*yy)
    rsd[t,1,1]<- 0.0
    rsd[t,2,1]<- -yy*zz*zz
    rsd[t,1,2]<- -yy*zz*zz
    rsd[t,2,2]<- 2.0*x[1]*yy*yy*zz*zz*zz
    rsd[t,3,1]<- t*x[2]*yy*zz*zz
    rsd[t,1,3]<- t*x[2]*yy*zz*zz
    rsd[t,3,2]<- t*x[1]*yy*zz*zz*(1-2*x[2]*yy*zz)
    rsd[t,2,3]<- t*x[1]*yy*zz*zz*(1-2*x[2]*yy*zz)
##    rsd[t,3,3]<- 2*t*t*x[1]*x[2]*x[2]*yy*yy*zz*zz*zz
    rsd[t,3,3]<- -t*t*x[1]*x[2]*yy*zz*zz*(1-2*yy*zz*x[2])
    return(rsd)
}


hobbs.h <- function(x) { ## compute Hessian
#   cat("Hessian not yet available\n")
#   return(NULL)
    H<-matrix(0,3,3)
    res<-hobbs.res(x)
    jj<-hobbs.jac(x)
    rsd<-hobbs.rsd(x)
##    H<-2.0*(t(res) %*% rsd + t(jj) %*% jj)
    for (j in 1:3) {
       for (k in 1:3) {
          for (i in 1:12) {
             H[j,k]<-H[j,k]+res[i]*rsd[i,j,k]
          }
       }
    }
    H<-2*(H + t(jj) %*% jj)
    return(H)
}


hobbsrsd.tst<-function(x) { # test rsd calculations
   hh<-1e-7 # use this for delta for derivatives
   Ja<-hobbs.jac(x)
   rsd<-hobbs.rsd(x)
   x1<-x+c(hh,0,0)
   x2<-x+c(0,hh,0)
   x3<-x+c(0,0,hh)
   Ja1<-hobbs.jac(x1)
   Ja2<-hobbs.jac(x2)
   Ja3<-hobbs.jac(x3)
   cat("w.r.t. x1 ")
   print(maxard((Ja1-Ja)/hh,rsd[,,1] ))
   cat("w.r.t. x2 ")
   print(maxard((Ja2-Ja)/hh,rsd[,,2] ))
   cat("w.r.t. x3 ")
   print(maxard((Ja3-Ja)/hh,rsd[,,3] ))
}


hobbs.doc <- function() { ## documentation for hobbs
   cat("One generalization of the Rosenbrock banana valley function (n parameters)\n")
   ## How should we do the documentation output?
}


hobbs.setup <- function(n=NULL, dotdat=NULL) {
  # if (is.null(gs) ) { gs<-100.0 } # set the scaling
  # if ( is.null(n) ) {
  #    n <- readline("Order of problem (n):")
  # }
   n<-3 # fixed for Hobbs, as is m=12
   x<-rep(2,n)
   lower<-rep(-100.0, n)
   upper<-rep(100.0, n)
   bdmsk<-rep(1,n)
   if (! is.null(dotdat) ) {
       fargs<-paste("gs=",gs,sep='') # ?? still need to do this nicely
   } else { fargs<-NULL }
   gsu<-list(x=x,lower=lower,upper=upper,bdmsk=bdmsk,fargs=fargs)
   return(gsu)

}

hobbs.fgh <- function(x) { # all 3 for trust method
         stopifnot(is.numeric(x))
         stopifnot(length(x) == 3)
         f<-hobbs.f(x)
         g<-hobbs.g(x)
         B<-hobbs.h(x)
         list(value = f, gradient = g, hessian = B)
}

hobbs.test<-function() { # test grad
require(numDeriv)
xvec<-rep(2,3)
cat("test hobbs with gs=100\n")
cat("function=",hobbs.f(xvec,),"\n")
gra<-hobbs.g(xvec,0)
cat("gra:")
print(gra)
gra0<-hobbs.g0(xvec,0)
cat("gra0:")
print(gra0)
grn<-grad(hobbs.f,xvec,0)
cat("grn:")
print(grn)
}

hobbs.time<-function() { # test grad
require(numDeriv)
 for (i in 1:6) {
  n<-2500*i
  xvec<-rep(2,n)
  cat("\n\n\n")
  cat("test hobbs with gs=100 and n=",n,"\n")
  cat("Function \n")
  print(  system.time(hobbs.f(xvec,)))
  cat("gra \n")
  print(system.time(hobbs.g(xvec,0)))
  cat("gra0\n")
  print(system.time(hobbs.g0(xvec,0)))
  cat("grn\n")
  print(system.time(grad(hobbs.f,xvec,0)))
 }
}

hobbs.examples<-function() { #test with different methods
   source("hobbs.R")
   x0<-c(1,1,1)
   ansbfgs0n<-optim(x0,hobbs.f,gr=NULL,method='BFGS')
   ansbfgs0a<-optim(x0,hobbs.f,hobbs.g,method='BFGS')
   xstar<-c(196.1770491, 49.0906349, 0.3135749)
   fstar<-hobbs.f(xstar)
   sc1<-c(100, 10, .1)
   xstars<-xstar/sc1
   xstars
   ansbfgs0ass<-optim(x0s,hobbs.f,hobbs.g,method='BFGS',control=list(parscale=sc1))
   ansbfgs0as<-optim(x0,hobbs.f,hobbs.g,method='BFGS',control=list(parscale=sc1))
   ansbfgs0ass<-optim(x0s,hobbs.f,hobbs.g,method='BFGS',control=list(parscale=sc1, trace=3))
   ansbfgs0n<-optim(x0,hobbs.f,gr=NULL,method='BFGS', control=list(trace=3))
   library(ucminf)
   library(BB)
   ansucminf0n<-ucminf(x0, hobbs.f,gr=NULL)
   ansucminf0n
   ansucminf0a<-ucminf(x0, hobbs.f,hobbs.g)
   ansucminf0a
   ansspg0n<-spg(x0,hobbs.f)
   ansspg0n
   ansspg0nx<-spg(x0,hobbs.f,control=list(maxit=50000))
   ansspg0nx
   ansspg0nxx<-spg(x0,hobbs.f,control=list(maxit=50000, maxfeval=1000000))
   ansspg0nxx
   anstrust<-trust(hobbs.fgh, x0, .1, 100, blather=TRUE)
   anstrust
   ansnlm0n<-nlm(hobbs.f,x0)
   ansnlm0n
   ansnlminb0<-nlminb(x0,hobbs.f)
   ansnlminb1<-nlminb(x0,hobbs.f, hobbs.g)
   ansnlminb2<-nlminb(x0,hobbs.f, hobbs.g,hobbs.h)
   ansnlminb0
   ansnlminb1
   ansnlminb2
   anspowell<-powell(x0,hobbs.f)
   anspowell
   ansocg0a<-optim(x0,hobbs.f,hobbs.g,method='CG')
   ansocg0a
   ansocg0n<-optim(x0,hobbs.f,method='CG')
   ansocg0n
   ansonm<-optim(x0,hobbs.f)
   ansonm
   ansolbfgsb0n<-optim(x0,hobbs.f,method='L-BFGS-B')
   ansolbfgsb0n
   ansosann0n<-optim(x0,hobbs.f,method='SANN')
   ansosann0n
   cat("nls tests\n")
   y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 50.156, 62.948,75.995, 91.972)
   t<-1:12
   xx<-c(200,10,.1)
   ansxx<-nls(y~x1/(1+x2*exp(-x3*t)), data=list(y,t), start=list(x1=xx[1],x2=xx[2],x3=xx[3]))
   xy<-c(1,1,1)
   ansxy<-nls(y~x1/(1+x2*exp(-x3*t)), data=list(y,t), start=list(x1=xy[1],x2=xy[2],x3=xy[3]))
}

require(snewton)
x0 <- c(200, 50, .3)
cat("Start for Hobbs:")
print(x0)
solx0 <- snewton(x0, hobbs.f, hobbs.g, hobbs.h)
print(solx0)
print(eigen(solx0$Hess)$values)

x1s <- c(100, 10, .1)
cat("Start for Hobbs:")
print(x1s)
solx1s <- snewton(x1s, hobbs.f, hobbs.g, hobbs.h, control=list(trace=2))
print(solx1s)
print(eigen(solx1s$Hess)$values)


x1 <- c(1, 1, 1)
cat("Start for Hobbs:")
print(x1)
solx1 <- snewton(x1, hobbs.f, hobbs.g, hobbs.h, control=list(trace=2))
print(solx1)
print(eigen(solx1$Hess)$values)

