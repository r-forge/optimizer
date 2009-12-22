###################################################
## badeval.R -- function and examples where we get inadmissible fn
rm(list=ls())
library(optimx)

if(!require("BB"))    stop("this requires package BB.")
if(!require(numDeriv))stop("this requires package numDeriv.")
if(!require("setRNG"))stop("this requires setRNG.")
if(!require("Rcgmin"))stop("this requires Rcgmin.")
if(!require("minqa"))stop("this requires minqa.")
if(!require("ucminf"))stop("this requires ucminf.")

# Replaced April 7, 2008, with setRNG to ensure rng and normal generators are set too.
test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=1236)
setRNG(test.rng)

## Functions for the tests **


badeval1.f <- function(x){
	n <- length(x)
	fval<-1
	for (i in 1:n) {
		fval=fval*sqrt(x[i]+i)
	}
	fval<-fval*fval
	return(fval)
}

##-----------------------------------------------------##
badeval1.cf <- function(x){
# This one has a check for badness
	n <- length(x)
	iseq<-seq(1:n)
	if (any((x[iseq]+iseq) < 0)) {
		fval<-0.5*.Machine$double.xmax
	} else {
		fval<-1
		for (i in 1:n) {
			fval=fval*sqrt(x[i]+i)
		}
		fval<-fval*fval
	}
	return(fval)
}
##-----------------------------------------------------##

badeval.f <- function(x){
	n <- length(x)
	iseq<-seq(1:n)
	fvec<-sqrt(x[iseq]+iseq)
	fval<-as.numeric(crossprod(fvec))
	return(fval)
}

##-----------------------------------------------------##
badeval.cf <- function(x){
# This one has a check for badness
	n <- length(x)
	iseq<-seq(1:n)
	if (any((x[iseq]+iseq) < 0)) {
		fval<-0.5*.Machine$double.xmax
	} else {
		fvec<-sqrt(x[iseq]+iseq)
		fval<-as.numeric(crossprod(fvec))
	}
	return(fval)
}
##-----------------------------------------------------##
 
##########
p0 <- abs(rnorm(2,sd=5)) ## start feasible
# ans.badevalc<-optimx(p0,fn=badeval.cf,control=list(all.methods=TRUE, save.failures=TRUE, trace=2))
ans.badevalc<-optimx(p0,fn=badeval.cf,control=list(all.methods=TRUE, save.failures=TRUE))
cat("\n")
cat("Problem badevalc -- check for inadmissible within function\n")
print(ans.badevalc)

tmpstr<-readline("Continue?")

iseq<-seq(1:2)
lower<--iseq
upper<-c(Inf,Inf)

# ans.badevalcb<-optimx(p0,fn=badeval.cf,gr=NULL, lower=lower, upper=upper, control=list(all.methods=TRUE, save.failures=TRUE, trace=2))
ans.badevalcb<-optimx(p0,fn=badeval.cf,gr=NULL, lower=lower, upper=upper, control=list(all.methods=TRUE, save.failures=TRUE))
cat("\n")
cat("Problem badevalcb - check for inadmissible within function AND bounds on parameters\n")
print(ans.badevalcb)

tmpstr<-readline("Continue?")

# ans.badeval<-optimx(p0,fn=badeval.f,control=list(all.methods=TRUE, save.failures=TRUE, trace=2))
ans.badeval<-optimx(p0,fn=badeval.f,control=list(all.methods=TRUE, save.failures=TRUE))
cat("\n")
cat("Problem badeval -- no checks\n")
print(ans.badeval)

