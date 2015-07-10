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

ufn <- function(x, parscale=c(1,1,1)){
    val <- hobbs.f(x*parscale)
}

cat("standard test using optim\n")
start <- c(300, 50, .3)
ps <- c(100, 10, .1)

anm <- optim(start, hobbs.f, control=list(trace=1))
anms <- optim(start, hobbs.f, control=list(trace=1, parscale=ps))

cat("no scaling\n")
print(anm)
cat("scaled using ps:")
print(ps)
print(anms)

tmp <- readline("Now try ufn scaling - first normal")

require(dfoptim)
anmk <- nmk(start, hobbs.f, control=list(trace=TRUE))
tmp <- readline("ufn unscaled")
anmku <- nmk(start, ufn, control=list(trace=TRUE), parscale=c(1,1,1))
starts <- start/ps
tmp <- readline("ufn scaling using starts")
anmks <- nmk(starts, ufn, control=list(trace=TRUE), parscale=ps)
cat("no scaling\n")
print(anmk)
cat("no scaling via c(1,1,1) in ps\n")
print(anmku)
cat("scaled using ps:")
print(ps)
print(anmks)


