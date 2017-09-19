# Test default line search  tdlsch.R
rm(list=ls())

wood.f <- function(x){
  res <- 100*(x[1]^2-x[2])^2+(1-x[1])^2+90*(x[3]^2-x[4])^2+(1-x[3])^2+
    10.1*((1-x[2])^2+(1-x[4])^2)+19.8*(1-x[2])*(1-x[4])
  return(res)
}
wood.g <- function(x){
  g1 <- 400*x[1]^3-400*x[1]*x[2]+2*x[1]-2
  g2 <- -200*x[1]^2+220.2*x[2]+19.8*x[4]-40
  g3 <- 360*x[3]^3-360*x[3]*x[4]+2*x[3]-2
  g4 <- -180*x[3]^2+200.2*x[4]+19.8*x[2]-40
  return(c(g1,g2,g3,g4))
}

#hessian:
wood.h <- function(x){
  h11 <- 1200*x[1]^2-400*x[2]+2;    h12 <- -400*x[1]; h13 <- h14 <- 0
  h22 <- 220.2; h23 <- 0;    h24 <- 19.8
  h33 <- 1080*x[3]^2-360*x[4]+2;    h34 <- -360*x[3]
  h44 <- 200.2
  H <- matrix(c(h11,h12,h13,h14,h12,h22,h23,h24,
                h13,h23,h33,h34,h14,h24,h34,h44),ncol=4)
  return(H)
}



rbk.f <- function(x){
  return(100*(x[2] - x[1]*x[1])^2 + (1-x[1])^2)
}
rbk.g <- function(x){
  return(c(-400*x[1]*(x[2] - x[1]*x[1]) - 2*(1-x[1]), 200*(x[2] - x[1]*x[1])))
}

#Hessian
rbk.h <- function(x) {
  a11 <- 2 - 400*x[2] + 1200*x[1]*x[1]; a21 <- -400*x[1]
  return(matrix(c(a11, a21, a21, 200), 2, 2))
}


control <- list(defstep=1,  trace = 0,
                maxit = 500,
                maxfevals = npar*500,
                epstol = .Machine$double.eps,
                stepdec = 0.2, 
                stepmax = 5,
                stepmin = 0,
                smult = 1.5,
                fmult = -0.25,
                offset = 100.0,
                bigval = .Machine$double.xmax*0.01stepmin=0) # these are interesting for Wood fn.

nf <- 0

lnsrch<-function(fn, fbest, xc, d, grv, ...) { # Line search using internal optimize()
  flsch<-function(st) {
    # computes the function value at stepsize st on line (xc + gm*d)
    # Essentially flsch(st)
    # gm: step size
    # fn: objective function
    # xc: base set of parameters
    # d : search direction
#    nf <- nf +1
    fval<-fn(xc+st*d,...)
    fval
  }

  lout<-fminsfp(flsch,0,control$defstep, ...)
  cat("lnsrch lout:")
  print(lout)
  rlout <- lout$xmin
#  cat("structure of rlout")
#  print(str(rlout))
  attr(rlout, "Fval") <- lout$fmin
  attr(rlout, "fcount") <- (lout$niter + 1) # fevals is iterations + 1
  rlout # Note: returns stepsize, not x
} # end default line search

cat("Rosenbrock:\n")

xr <- c(-1.2, 1)
dr <- solve(rbk.h(xr), -rbk.g(xr))

f0 <- rbk.f(xr)+1
cat("f0=",f0,"\n")
cat("Initil search direction:")
print(dr)


tr <- lnsrch(rbk.f, f0, xr, dr)
cat("Solution from default line search:")
print(tr)

ff<-rep(NA,100)

xbits <- seq(1:100)/60

for (i in 1:100){ 
  y <- xr+xbits[i]*dr
  ff[i]<- rbk.f(y)
}
plot(xbits, ff)
cat("minimum position from search on plot points:")
print(xbits[which(ff == min(ff))])

tmp <- readline("Continue?")

cat("Wood Function:\n")
xw <- c(-3,-1,-3,-1) # Wood standard start

dw <- solve(wood.h(xw), -wood.g(xw))

f0 <- wood.f(xw)+1
cat("f0=",f0,"\n")
cat("Initil search direction:")
print(dw)


tw <- lnsrch(wood.f, f0, xw, dw)
cat("Solution from default line search:")
print(tw)

ffw<-rep(NA,100)

xwbits <- seq(1:100)/60

for (i in 1:100){ 
  y <- xw+xwbits[i]*dw
  ffw[i]<- wood.f(y)
}
plot(xwbits, ffw)
cat("minimum position from search on plot points:")
print(xwbits[which(ffw == min(ffw))])

