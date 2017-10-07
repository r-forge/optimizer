## cyq.R -- Fletcher's Chebyquad problem
## @knitr cyq.prb
# This is file cyq.prb
probname <- "cyq"
probdesc <- "
 Ref: Fletcher, R. (1965) Function minimization without calculating 
 derivatives -- a review, Computer J., 8, 33-41.
 More et al problem 35

Note that this has m >= n cases, but often uses m=n in tests.

"

#- Note: environment / list "pe" must already exist
if (exists("pe")) { 
      rm("pe")  
  }

pe <- new.env()
pe$kf <- 0
pe$kg <- 0
pe$kjac <- 0
pe$kres <- 0

#- nls format expression
## cyq.formula <- ??

#- setup
# cyq.setup<-function() {

# Note we do not have all components here e.g., .jsd, .h

cyq.f <- function (x) {
  rv<-cyq.res(x)
  f<-sum(rv*rv)
}

#- gradient
cyq.g <- function (x) {
   cj<-cyq.jac(x)
   rv<-cyq.res(x)
   gg<- as.vector(2.0* rv %*% cj)
}

#- hessian
cyq.h <- function(x) {
#- THIS IS NOT COMPLETE??
  stop("cyq.h IS NOT COMPLETE")
  res<-cyq.res(x)
  JJ<-cyq.jac(x)
  pe$khess <- pe$khess + 1
  H <- t(JJ) %*% JJ
}


#- residual
cyq.res <- function (x) {
# Fletcher's chebyquad function m = n -- residuals 
   n<-length(x)
   res<-rep(0,n) # initialize
   for (i in 1:n) { #loop over resids
     rr<-0.0
     for (k in 1:n) {
	z7<-1.0
	z2<-2.0*x[k]-1.0
        z8<-z2
        j<-1
        while (j<i) {
            z6<-z7
            z7<-z8
            z8<-2*z2*z7-z6 # recurrence to compute Chebyshev polynomial
            j<-j+1
        } # end recurrence loop
        rr<-rr+z8
      } # end loop on k
      rr<-rr/n
      if (2*trunc(i/2) == i) { rr <- rr + 1.0/(i*i - 1) }
      res[i]<-rr
    } # end loop on i
    res
}

#- Jacobian
cyq.jac<- function (x) {
#  Chebyquad Jacobian matrix
   n<-length(x)
   cj<-matrix(0.0, n, n)
   for (i in 1:n) { # loop over rows
     for (k in 1:n) { # loop over columns (parameters)
       z5<-0.0
       cj[i,k]<-2.0
       z8<-2.0*x[k]-1.0 
       z2<-z8
       z7<-1.0
       j<- 1
       while (j<i) { # recurrence loop
         z4<-z5
         z5<-cj[i,k]
         cj[i,k]<-4.0*z8+2.0*z2*z5-z4
         z6<-z7
         z7<-z8
         z8<-2.0*z2*z7-z6
         j<- j+1
       } # end recurrence loop
       cj[i,k]<-cj[i,k]/n
     } # end loop on k
   } # end loop on i
   cj
}

#- nleq -- ?? section if appropriate

#- test
## put example calls of the function, possibly including calls to 
# optimizations and nonlinear least squares etc.

cyq.setup <- function(n = NULL) { 
  cat("Fletcher chebyquad function in file cyq.R\n")
  if (is.null(n)) {
     n <- as.numeric(readline("order of problem (n) ="))
  }
  lower<-rep(-10.0, n)
  upper<-rep(10.0, n) 
  bdmsk<-rep(1, n) # free all parameters
    x<-1:n
    x<-x/(n+1.0) # Initial value suggested by Fletcher
  result<-list(x=x, lower=lower, upper=upper, bdmsk=bdmsk)
}

## Examples 

for (n in 2:9) { 
  cs <- cyq.setup(n)
  library(optimr)
  cat("x0=")
  x0 <- cs$x
  print(x0)
  cat("\n")
  cat("Chebyquad function at x0 = ",cyq.f(x0),"\n")  
  opmcyq <- opm(x0, cyq.f, cyq.g, method="ALL")
  print(summary(opmcyq, order=value, par.select=1:4))
}


#- End cyq.prb   
