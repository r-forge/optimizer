# Test backtrack line search  tdlsch.R
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
xw <- c(-3,-1,-3,-1) # Wood standard start

dw <- solve(wood.h(xw), -wood.g(xw))




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


xr <- c(-1.2, 1)
dr <- solve(rbk.h(xr), -rbk.g(xr))

cat("Initil search direction:")
print(dr)

control <- list(stepdec=0.2, offset=100., trace=2, acctol=0.0001)
nf <- 0
grv <- rbk.g(xr)
cat("Initial Rbk fn=",rbk.f(xr),"\n")
cat("Initial gradient =")
print(grv)

lnsrch<-function(fn, fbest, xc, d, grv, ...) { # backtrack line search
    flsch<-function(st) {
      # computes the function value at stepsize st on line (xc + gm*d)
      # Essentially flsch(st)
      # gm: step size
      # fn: objective function
      # xc: base set of parameters
      # d : search direction
      nf <- nf +1
      fval<-fn(xc+st*d,...)
      fval
    }
    st <- 1.0
    gproj <- as.numeric(crossprod(grv, d) )
    cat("gradproj =",gproj,"\n")
    repeat {
      xnew <- xc + st*d # new point
      if ((control$offset+st) == (control$offset)) { # no better parameters
        st <- 0
        rlout <- st
        attr(rlout,"Fval")<-fbest # Assume we pass this in
        return(rlout)
      }
      fval <- flsch(xnew, ...)
      if (control$trace > 1) cat("Step = ",st," fval = ",fval,"\n")
      if (fval <= fbest + control$acctol*st*gproj) break
      st <- control$stepdec*st # new step
    }
    rlout <- st
    attr(rlout, "Fval")<- fval
    rlout
  } # end backtrack line search
  
f0 <- rbk.f(xr)+1
f0
cat("f0=",f0,"\n")

tr <- lnsrch(rbk.f, f0, xr, dr, grv)
print(tr)

