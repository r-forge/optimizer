rm(list=ls())
require(snewton)
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


cat("This appears to fail\n")
x1 <- c(1, 1, 1)
cat("Start for Hobbs:")
print(x1)
ftest <- try(solx1 <- snewton(x1, hobbs.f, hobbs.g, hobbs.h, control=list(trace=2)))
if (class(ftest) != "try-error") {
  print(solx1)
  cat("Eigenvalues of returned Hessian: ")
  print(eigen(solx1$Hess)$values)
}
ftestm <- try(solx1m <- snewtonm(x1, hobbs.f, hobbs.g, hobbs.h, control=list(trace=2)))
if (class(ftestm) != "try-error") {
  print(solx1m)
  cat("Eigenvalues of returned Hessian: ")
  print(eigen(solx1m$Hess)$values)
}


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

cat("This appears to fail\n")
x1 <- c(1, 1, 1)
cat("Start for Hobbs:")
print(x1)
ftest <- try(solx1 <- snewton(x1, hobbs.f, hobbs.g, hobbs.h, control=list(trace=2)))
if (class(ftest) != "try-error") {
   print(solx1)
   print(eigen(solx1$Hess)$values)
}
# we can also use nlm and nlminb
#??

xtemp <- c(12.156272,  2.321825,  4.166667)
fxtemp <- hobbs.f(xtemp)
gxtemp <- hobbs.g(xtemp)
fxtemp
gxtemp


xb <- c(12.156272,  2.321825,  4.166667)
fxb <- hobbs.f(xb)
fxb

xc <- c(11.716841,  2.317431,  4.022304)
print(hobbs.f(xc))
