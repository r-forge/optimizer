## @knitr BroydenTestR.prb
probname <- "BroydenTestR"
probdesc <- "This problem found in package BB. It appears to be similar, 
and may compute to the same value as, either Broyden banded or Broyden
tridiagonal problems. refs.??
"
probtype <- "nls" # nls offers more possibilities than function minimization
#-?? should we list all types available in THIS file?
#- Other types??  "fmin", "ssmin" (for resid ss)
mformula <- NA
#- The above is a "formula". But we could also possibl also use a character string or expression,
#- but need to know how to do so carefully. ??
##- ?? mdata=BroydenTestR # This would work if we have a data package, but not useful independently
#- This is the local data frame. We need to rename this when we put it into a global structure
#- Now add various functions for optimization or nls-by-function
#- ?? We could have problems with quotation marks WITHIN the code.

BroydenTestR.f <- function(x) {
   n <- length(x)
   res <- rep(NA, n)
   res[1] <- ((3 - 0.5*x[1]) * x[1]) - 2*x[2] + 1
   tnm1 <- 2:(n-1)
   res[tnm1] <- ((3 - 0.5*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
   res[n] <- ((3 - 0.5*x[n]) * x[n]) - x[n-1] + 1
   sum(res*res)
}


#- resfn
BroydenTestR.res <- function(x) {
  res <- rep(NA, n)
  res[1] <- ((3 - 0.5*x[1]) * x[1]) - 2*x[2] + 1
   tnm1 <- 2:(n-1)
   res[tnm1] <- ((3 - 0.5*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
   res[n] <- ((3 - 0.5*x[n]) * x[n]) - x[n-1] + 1
   return(res)
}

#- jacfn
# BroydenTestR - Jacobian
BroydenTestR.jac <- function(x) {
   n <- length(x)

   # ?? to be fixed
   jj<-matrix(NA,n,n) # gradient set to NA to start with
   f <- rep(NA, n)
   f[1] <- ((3 - 0.5*x[1]) * x[1]) - 2*x[2] + 1
   tnm1 <- 2:(n-1)
   f[tnm1] <- ((3 - 0.5*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
   f[n] <- ((3 - 0.5*x[n]) * x[n]) - x[n-1] + 1
   sum(f*f)
   return(NA)
}

#- hessfn
BroydenTestR.h <- function(x) {
## JJ<-BroydenTestR.jac(x)
## H <- t(JJ) %*% JJ
## res<-BroydenTestR.res(x)
}

#- gr
BroydenTestR.g<-function(x) {
   n <- length(x)
   gg<-rep(NA,n) # gradient set to NA to start with
#   res <- rep(NA, n)
#   res[1] <- ((3 - 0.5*x[1]) * x[1]) - 2*x[2] + 1
#   res[tnm1] <- ((3 - 0.5*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
#   res[n] <- ((3 - 0.5*x[n]) * x[n]) - x[n-1] + 1
   ## Need to get partial derivatives of res w.r.t each x[j] ??
   # ??
  gg[1]<- -2 + 2*x[1] + 4*x[3] + (6 - 2.0*x[1])*(1 - 2*x[2] + x[1]*(3 - 0.5*x[1])) - 2*x[2]*(3 - 0.5*x[2])
  gg[2] <- -6 + 4*x[4] + 10*x[2] + (6 - 2.0*x[2])*(1 - x[1] - 2*x[3] + x[2]*(3 - 0.5*x[2])) - 4*x[1]*(3 - 0.5*x[1]) - 2*x[3]*(3 - 0.5*x[3])
 tnm2 <- 3:(n-2)
 gg[tnm2]<- -6 + 4*x[tnm2-2] + 4*x[tnm2+2] + 10*x[tnm2] + (6 - 2.0*x[tnm2])*(1 - x[tnm2-1] - 2*x[tnm2+1] 
             + x[tnm2]*(3 - 0.5*x[tnm2])) - 4*x[tnm2-1]*(3 - 0.5*x[tnm2-1]) - 2*x[tnm2+1]*(3 - 0.5*x[tnm2+1])
 gg[n-1] <- -6 + 4*x[n-3] + 10*x[n-1] + (6 - 2.0*x[n-1])*(1 - x[n-2] - 2*x[n] + x[n-1]*(3 - 0.5*x[n-1])) - 4*x[n-2]*(3 - 0.5*x[n-2]) - 2*x[n]*(3 - 0.5*x[n])

##          -6 + 4*x[2] +   10*x[4] +   (6 - 2.0*x[4])*(1 - x[3] - 2*x[5] + x[4]*(3 - 0.500000000000000*x[4])) - 4*x[3]*(3 - 0.500000000000000*x[3]) - 2*x[5]*(3 - 0.500000000000000*x[5])

 gg[n] <- -4 + 4*x[n-2] + 8*x[n] + (6 - 2.0*x[n])*(1 - x[n-1] + x[n]*(3 - 0.5*x[n])) - 4*x[n-1]*(3 - 0.5*x[n-1])
 return(gg)
}

## c( -2 + 2*x[1] + 4*x[3] + (6 - 2.00000000000000*x[1])*(1 - 2*x[2] + x[1]*(3 - 0.500000000000000*x[1])) - 2*x[2]*(3 - 0.500000000000000*x[2]),-6 + 4*x[4] + 10*x[2] + (6 - 2.00000000000000*x[2])*(1 - x[1] - 2*x[3] + x[2]*(3 - 0.500000000000000*x[2])) - 4*x[1]*(3 - 0.500000000000000*x[1]) - 2*x[3]*(3 - 0.500000000000000*x[3]),-6 + 4*x[1] + 4*x[5] + 10*x[3] + (6 - 2.00000000000000*x[3])*(1 - x[2] - 2*x[4] + x[3]*(3 - 0.500000000000000*x[3])) - 4*x[2]*(3 - 0.500000000000000*x[2]) - 2*x[4]*(3 - 0.500000000000000*x[4]),-6 + 4*x[2] + 10*x[4] + (6 - 2.00000000000000*x[4])*(1 - x[3] - 2*x[5] + x[4]*(3 - 0.500000000000000*x[4])) - 4*x[3]*(3 - 0.500000000000000*x[3]) - 2*x[5]*(3 - 0.500000000000000*x[5]),-4 + 4*x[3] + 8*x[5] + (6 - 2.00000000000000*x[5])*(1 - x[4] + x[5]*(3 - 0.500000000000000*x[5])) - 4*x[4]*(3 - 0.500000000000000*x[4]) )

#- setup
#- BroydenTestR.setup<-function() {
BroydenTestR.df <- NA
# out<-list(starts=starts, df=BroydenTestR.df)
#- return(out)
#- }

#- broydt.setup <- function(n=NULL, choice=NULL) {
#- ?? need to name the parameters and set up starts array
   cat("broydt.setup -- str of choice:")
   str(choice)
   if (is.null(choice)) choice<-1
   if(choice == 1) {
	start<-rep(-1,n)
   } else {
	if (choice == 2) {
		start<-rep(0,n)
	} else {
		warning("Inadmissible choice in broydt.setup")
		start<-NULL
	}
   } # end choice setup
options(digits=12)
if(!require("BB"))stop("this test requires package BB.")
if(!require("setRNG"))stop("this test requires setRNG.")

# Use a preset seed so test values are reproducable. 
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
old.seed <- setRNG(test.rng)

##########
cat("BB test broydt.f ...\n")

broydt.f <- function(x) {
n <- length(x)
f <- rep(NA, n)
f[1] <- ((3 - 0.5*x[1]) * x[1]) - 2*x[2] + 1
tnm1 <- 2:(n-1)
f[tnm1] <- ((3 - 0.5*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
f[n] <- ((3 - 0.5*x[n]) * x[n]) - x[n-1] + 1
sum(f*f)
}

p0 <- rnorm(50, sd=1)
system.time(ans.spg <- spg(par=p0, fn=broydt.f, control=list(maxit=10000)))[1]
 
z <- sum(ans.spg$par)
# good   <-   13.6791297530 converged ?
good   <-   3.393768541175171
print(z, digits=16)
if(any(abs(good - z) > 5e-1)) stop("BB test broydt.f a FAILED")
### =========================================================================== 


#- broydt.nlm <- function(x) {
#-     broydtnlm<-broydt.f(x)
#-     attr(broydtnlm, "gradient") <- broydt.g(x)
#-     attr(broydtnlm, "hessian") <- NULL
#-     return(broydtnlm)
#- }





BroydenTestR.test<-function() {
}   

