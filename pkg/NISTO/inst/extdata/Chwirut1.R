# Optimization test function chwirut1
# chwirut1 from NISTnls
# ??ref...


chwirut1.f <- function(x) {
   res<-chwirut1.res(x)
   f<-sum(res*res)
}

chwirut1.res <- function(x) {
   xx<-Chwirut1$x # note caps!
   yy<-Chwirut1$y
   res <- rep(NA, length(xx))
   b1<-x[1]
   b2<-x[2]
   b3<-x[3]
   res<- exp(-b1*xx)/(b2+b3*xx) - yy
   return(res)
}

# Chwirut1 function - Jacobian
chwirut1.jac <- function(x) {
## stop("not defined")
xx<-Chwirut1$x # note caps!
J <- matrix(0, length(xx), 3) # define the size of the Jacobian
   b1<-x[1]
   b2<-x[2]
   b3<-x[3]
#   res<-b1*(b2+xx)**(-1/b3) - yy
   exx<-exp(-b1*xx)
   den<-(b2+b3*xx)
   expr1 <- b2 + xx
    expr3 <- -1/b3
    expr4 <- expr1^expr3
    J[ , 1] <- -xx*exx/den
    J[ , 2] <- -exx/(den*den)
    J[ , 3] <- -exx*xx/(den*den)

return(J)
}

chwirut1.h <- function(x) {
   JJ<-chwirut1.jac(x)
   H <- t(JJ) %*% JJ
   res<-chwirut1.res(x)
stop("not defined")

}

chwirut1.g<-function(x) {
   JJ<-chwirut1.jac(x)
   res<-chwirut1.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

chwirut1.fgh<-function(x) {
   f<-chwirut1.f(x)
   g<-chwirut1.g(x)
   H<-chwirut1.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

chwirut1.setup<-function() {
   library(NISTnls) # get parent collection
   data(Chwirut1) # and load up the data into x and y
   start<-c(0.1, 0.01, 0.02)
}

chwirut1.test<-function() {
# Here is where we put the test code, but we need to save stuff off for later use
#    especially if there are multiple outputs.
# Better to call optimx, but need to capture output (sink()?)

 start1 = c( 0.1,  0.01,  0.02)
 ansnm1<-optim(start1, chwirut1.f, control=list(trace=3, maxit=1000000))
 start2 = c( 0.15,  0.008,  0.010)
 ansnm2<-optim(start2, chwirut1.f, control=list(trace=3, maxit=1000000))
 startj<-c(1,1,1)
 ansnmj<-optim(startj, chwirut1.f, control=list(trace=3, maxit=1000000))
  
 start1 = c( 0.1,  0.01,  0.02)
 ansbfgsn1<-optim(start1, chwirut1.f, method='BFGS', control=list(trace=3, maxit=1000000))
 start2 = c( 0.15,  0.008,  0.010)
 ansbfgsn2<-optim(start2, chwirut1.f, method='BFGS', control=list(trace=3, maxit=1000000))
 startj<-c(1,1,1)
 ansbfgsnj<-optim(startj, chwirut1.f, method='BFGS', control=list(trace=3, maxit=1000000))

 start1 = c( 0.1,  0.01,  0.02)
 anscgfrn1<-optim(start1, chwirut1.f, method='CG', control=list(trace=3, maxit=1000000))
 start2 = c( 0.15,  0.008,  0.010)
 anscgfrn2<-optim(start2, chwirut1.f, method='CG', control=list(trace=3, maxit=1000000))
 startj<-c(1,1,1)
 anscgfrnj<-optim(startj, chwirut1.f, method='CG', control=list(trace=3, maxit=1000000))
  
  
}
