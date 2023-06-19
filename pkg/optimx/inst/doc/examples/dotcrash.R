# dotcrash.R -- dotargs have a name that clashes with internals 
# In particular numDeriv grad has "x" as an argument.
# Might be that pracma works, since x0 is used instead of x??
rm(list=ls())
sqmod<-function(z, x){
   nn<-length(z)
   yy<-x^(1:nn)
   f<-sum((yy-z)^2)
#   cat("Fv=",f," at ")
#   print(z)
   f
}
sqmod.g <- function(z, x){
   nn<-length(z)
   yy<-x^(1:nn)
   gg<- 2*(z - yy)
}

sessionInfo()
nn<-2
st<-rep(0.5, nn)
ginp <- function(fn, vv, ...){
  gg <- pracma::grad(fn, x0=vv, ...)
  gg
}
ginn <- function(fn, vv, ...){
  gg <- numDeriv::grad(fn, vv, ...)
  gg
}
cat("analytic:")
print(sqmod.g(st, x=2))
cat("by pracma::grad:")
print(ginp(sqmod, st, x=2))
cat("by numDeriv::grad:")
print(ginn(sqmod, st, x=2))

# Workarounds?
sqmod1 <- function(z){ sqmod(z, x=x) }
dots <- list(x=2)
sqmod2 <- function(z){ sqmod(z, unlist(dots)) }

