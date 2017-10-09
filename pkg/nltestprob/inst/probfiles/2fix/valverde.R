## @knitr ##Valverde-Dani.prb
# This is file ##Valverde-Dani.prb
rm(list=ls())
probname <- "##Valverde-Dani"
probdesc <- "Put your description in double quotes.
"

#- Note: environment / list "counters" must already exist

if (exists("pe")) { 
      rm("pe")  
  }

pe <- new.env()
pe$kf <- 0
pe$kg <- 0
pe$kjac <- 0
pe$kres <- 0

#- nls format expression
##Valverde.formula <- ( y ~ b1*x**b2 )

#- setup

## library("NISTnls", character.only=TRUE)
## mypdata <- eval(parse(text=data("Valverde-Dani")))
## Optimization test function dani
## ?? refs (put in .doc??)
## Nash and Walker-Smith (1987, 1989) ...


dani.f<- function(x, y=y, t=t){ # # dani weeds problem -- function
#    cat("dani.f, x:")
#    print(x)
#    print(t)
#    print(y)
#    cat("(abs(max(t))*x[3]):",(abs(max(t))*x[3]),"\n")
    if ((abs(max(t))*x[3]) > 50) { # check computability
       fbad<-.Machine$double.xmax
       return(fbad)
    }
    res<-dani.res(x, y=y, t=t)
    f<-sum(res*res)
}
dani.f2<- function(x, y=y, t=t){ # # dani weeds problem -- function
#    cat("dani.f, x:")
#    print(x)
#    print(t)
#    print(y)
#    cat("(abs(max(t))*x[3]):",(abs(max(t))*x[3]),"\n")
    res<-dani.res2(x, y=y, t=t)
    f<-sum(res*res)
}


dani.res<-function(x, y=y, t=t){ # dani weeds problem -- residual
# This variant uses looping
    m<-length(y)
#    cat(m,"  elements\n")
    if(length(x) != 3) stop("dani.res -- parameter vector n!=3")
    if((abs(max(t))*x[3]) > 50) {
       res<-rep(Inf,m)
    } else {
       res<-x[1]/(1+x[2]*exp(-x[3]*t)) - y
    }
}

dani.res2<-function(x, y=y, t=t){ # dani weeds problem -- residual
# This variant uses looping
    m<-length(y)
#    cat(m,"  elements\n")
    if(length(x) != 3) stop("dani.res -- parameter vector n!=3")
#    if((abs(max(t))*x[3]) > 50) {
#       res<-rep(Inf,m)
#    } else {
       res<-x[1]/(1+exp((x[2]-t)/x[3])) - y
#    }
}

dani.jac<-function(x, y=y, t=t){ # Jacobian of dani weeds problem
   m<-length(y)
   jj<-matrix(0.0, m, 3)
    yy<-exp(-x[3]*t)
    zz<-1.0/(1+x[2]*yy)
    ii<-seq(1,m)
     jj[ii,1] <- zz
     jj[ii,2] <- -x[1]*zz*zz*yy
     jj[ii,3] <- x[1]*zz*zz*yy*x[2]*t
   return(jj)
}

dani.g<-function(x, y=y, t=t){ # gradient of dani weeds problem
    # NOT EFFICIENT TO CALL AGAIN
    jj<-dani.jac(x,y=y,t=t)
    res<-dani.res(x,y=y,t=t)
    gg<-as.vector(2.*t(jj) %*% res)
    return(gg)
}


dani.rsd<-function(x, y=y, t=t) { # Jacobian second derivative
    rsd<-array(0.0, c(lenght(y),3,3))
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


dani.h <- function(x, y=y, t=t) { ## compute Hessian
#   cat("Hessian not yet available\n")
#   return(NULL)
    H<-matrix(0,3,3)
    res<-dani.res(x, y=y, t=t)
    jj<-dani.jac(x, y=y, t=t)
    rsd<-dani.rsd(x, y=y, t=t)
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


danirsd.tst<-function(x) { # test rsd calculations
   hh<-1e-7 # use this for delta for derivatives
   Ja<-dani.jac(x, y=y, t=t)
   rsd<-dani.rsd(x, y=y, t=t)
   x1<-x+c(hh,0,0)
   x2<-x+c(0,hh,0)
   x3<-x+c(0,0,hh)
   Ja1<-dani.jac(x1)
   Ja2<-dani.jac(x2)
   Ja3<-dani.jac(x3)
   cat("w.r.t. x1 ")
   print(maxard((Ja1-Ja)/hh,rsd[,,1] ))
   cat("w.r.t. x2 ")
   print(maxard((Ja2-Ja)/hh,rsd[,,2] ))
   cat("w.r.t. x3 ")
   print(maxard((Ja3-Ja)/hh,rsd[,,3] ))
}


dani.doc <- function() { ## documentation for dani
   cat("One generalization of the Rosenbrock banana valley function (n parameters)\n")
   ## How should we do the documentation output?
}


dani.setup <- function(n=NULL, dotdat=NULL) {
  # if (is.null(gs) ) { gs<-100.0 } # set the scaling
  # if ( is.null(n) ) {
  #    n <- readline("Order of problem (n):")
  # }
   n<-3 # fixed for dani, as is m=12
   y<-c(0.2000000,  0.1958350, 0.2914560, 0.6763628, 0.8494534, 0.9874526, 1.0477692)
   t<- c(0, 10,20, 40, 60, 90, 120)
   x<-rep(2,n)
   lower<-rep(-100.0, n)
   upper<-rep(100.0, n)
   bdmsk<-rep(1,n)
   if (! is.null(dotdat) ) {
       fargs<-paste("gs=",gs,sep='') # ?? still need to do this nicely
   } else { fargs<-NULL }
   gsu<-list(x=x,y=y, t=t,lower=lower,upper=upper,bdmsk=bdmsk,fargs=fargs)
   return(gsu)

}

dani.fgh <- function(x) { # all 3 for trust method
         stopifnot(is.numeric(x))
         stopifnot(length(x) == 3)
         f<-dani.f(x)
         g<-dani.g(x)
         B<-dani.h(x)
         list(value = f, gradient = g, hessian = B)
}

dani.test<-function() { # test grad
require(numDeriv)
xvec<-rep(2,3)
cat("test dani with gs=100\n")
cat("function=",dani.f(xvec,),"\n")
gra<-dani.g(xvec,0)
cat("gra:")
print(gra)
gra0<-dani.g0(xvec,0)
cat("gra0:")
print(gra0)
grn<-grad(dani.f,xvec,0)
cat("grn:")
print(grn)
}

dani.time<-function() { # test grad
require(numDeriv)
 for (i in 1:6) {
  n<-2500*i
  xvec<-rep(2,n)
  cat("\n\n\n")
  cat("test dani with gs=100 and n=",n,"\n")
  cat("Function \n")
  print(  system.time(dani.f(xvec,)))
  cat("gra \n")
  print(system.time(dani.g(xvec,0)))
  cat("gra0\n")
  print(system.time(dani.g0(xvec,0)))
  cat("grn\n")
  print(system.time(grad(dani.f,xvec,0)))
 }
}

source("dani.R")
   x0<-c(1,1,1)
   ansbfgs0n<-optim(x0,dani.f,gr=NULL,method='BFGS')
   ansbfgs0a<-optim(x0,dani.f,dani.g,method='BFGS')
   xstar<-c(196.1770491, 49.0906349, 0.3135749)
   fstar<-dani.f(xstar)
   sc1<-c(100, 10, .1)
   xstars<-xstar/sc1
   xstars
   ansbfgs0ass<-optim(x0s,dani.f,dani.g,method='BFGS',control=list(parscale=sc1))
   ansbfgs0as<-optim(x0,dani.f,dani.g,method='BFGS',control=list(parscale=sc1))
   ansbfgs0ass<-optim(x0s,dani.f,dani.g,method='BFGS',control=list(parscale=sc1, trace=3))
   ansbfgs0n<-optim(x0,dani.f,gr=NULL,method='BFGS', control=list(trace=3))
   library(ucminf)
   library(BB)
   ansucminf0n<-ucminf(x0, dani.f,gr=NULL)
   ansucminf0n
   ansucminf0a<-ucminf(x0, dani.f,dani.g)
   ansucminf0a
   ansspg0n<-spg(x0,dani.f)
   ansspg0n
   ansspg0nx<-spg(x0,dani.f,control=list(maxit=50000))
   ansspg0nx
   ansspg0nxx<-spg(x0,dani.f,control=list(maxit=50000, maxfeval=1000000))
   ansspg0nxx
   anstrust<-trust(dani.fgh, x0, .1, 100, blather=TRUE)
   anstrust
   ansnlm0n<-nlm(dani.f,x0)
   ansnlm0n
   ansnlminb0<-nlminb(x0,dani.f)
   ansnlminb1<-nlminb(x0,dani.f, dani.g)
   ansnlminb2<-nlminb(x0,dani.f, dani.g,dani.h)
   ansnlminb0
   ansnlminb1
   ansnlminb2
   anspowell<-powell(x0,dani.f)
   anspowell
   ansocg0a<-optim(x0,dani.f,dani.g,method='CG')
   ansocg0a
   ansocg0n<-optim(x0,dani.f,method='CG')
   ansocg0n
   ansonm<-optim(x0,dani.f)
   ansonm
   ansolbfgsb0n<-optim(x0,dani.f,method='L-BFGS-B')
   ansolbfgsb0n
   ansosann0n<-optim(x0,dani.f,method='SANN')
   ansosann0n
   cat("nls tests\n")
   y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 50.156, 62.948,75.995, 91.972)
   t<-1:12
   xx<-c(200,10,.1)
   ansxx<-nls(y~x1/(1+x2*exp(-x3*t)), data=list(y,t), start=list(x1=xx[1],x2=xx[2],x3=xx[3]))
   xy<-c(1,1,1)
   ansxy<-nls(y~x1/(1+x2*exp(-x3*t)), data=list(y,t), start=list(x1=xy[1],x2=xy[2],x3=xy[3]))






AMP.nls <- nls(AMP~SSlogis(Time,Asym, xmid, scal), data = concentrations,model=T)

AMP.nls <- nls(y~SSlogis(t,Asym, xmid, scal), data = concentrations,model=T)


