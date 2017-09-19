rm(list=ls())

  efn <- function(spar, ...) {
      # rely on pscale being defined in this enclosing environment
      par <- spar*pscale
      val <- fn(par, ...)
  }
  egr <- function(spar, ...) {
      par <- spar*pscale
      result <- gr(par, ...) * pscale
  }

  nlmfn <- function(spar, ...){
     f <- efn(spar, ...)
     g <- egr(spar, ...)
     attr(f,"gradient") <- g
     attr(f,"hessian") <- NULL # ?? maybe change later
     f
  }

fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
    x1 <- x[1]
    x2 <- x[2]
    c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
       200 *      (x2 - x1 * x1))
}

start <- c(-1.2,1)
fnscale <- 1
pscale <- rep(1,2)

fn <- fr
gr <- grr

print(efn(start*pscale))
print(egr(start*pscale))
print(nlmfn(start*pscale))

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
    y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
         38.558, 50.156, 62.948, 75.995, 91.972)
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


start <- c(300, 50, .3)
ps <- c(100, 10, .1)

pscale <- rep(1,3)
fn <- hobbs.f
gr <- hobbs.g

print(efn(start*pscale))
print(egr(start*pscale))
print(nlmfn(start*pscale))

sts <- start/ps

pscale <- ps
print(efn(sts))
print(egr(sts))
print(nlmfn(sts))
# check
require(numDeriv)
print(grad(efn, sts))

