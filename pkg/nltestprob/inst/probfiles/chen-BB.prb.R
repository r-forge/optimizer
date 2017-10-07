## @knitr ##Chen-BB.prb
# This is file ##Chen-BB.prb
probname <- "##Chen-BB"
probdesc <- "
This problem appears as one of the tests of the BB package for R on CRAN.

So far I (JN) have not found a source for this problem as at 20171007.
It appears to be a small-residual or nonlinear equations problem with
variable number of parameters.


The problem has only a functional form. I added the gradient.
It appears that there may be multiple solutions or near-solutions.

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
##Chen-BB.formula <- ( y ~ b1*x**b2 )
## Chen function ?? ref?

chen.f <- function(x) {
   v <- log(x) + exp(x)
   res <- (v - sqrt(v^2 + 5e-04))/2
   sum (res * res)
}

chen.g <- function(x) {
   v <- log(x) + exp(x)
   res <- (v - sqrt(v^2 + 5e-04))/2
   jj<-chen.jac(x)
   gg<- 2.0 * jj %*% as.vector(res)
   return(gg)
}

chen.res <- function(x) {
   v <- log(x) + exp(x)
   res <- (v - sqrt(v^2 + 5e-04))/2
}

chen.jac <- function(x) {
   n<-length(x)
   v <- log(x) + exp(x)
   jj<-matrix(0.0, n, n)
   for (i in 1:n) {
     jj[i,i] <- 0.5 * (1.0/x[i] + exp(x[i])) * (1.0-v[i]/sqrt(v[i]^2 + 5e-04))
   } 
   return(jj)
}




chen.setup <- function(n = 2) { 
  require(setRNG)
  test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
  old.seed <- setRNG(test.rng)
  st <- rexp(n)   
  lower<-rep(0.0, n) # Because log cannot take negative argument?
  upper<-rep(100.0, n)
  bdmsk<-rep(1,n)
  gsu<-list(st=st,lower=lower,upper=upper,bdmsk=bdmsk)
  gsu
}

## Examples

n <- 2
cset <- chen.setup(n)
cat("start =")
print(cset$st)
cBB2opm <- opm(cset$st, chen.f, chen.g, method="ALL",  lower=cset$lower, upper=cset$upper, control=list(trace=1))
summary(cBB2opm, order=value)




