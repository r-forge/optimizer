# Optimization test function chwirut2
# chwirut2 from NISTnls
# ??ref...


chwirut2.f <- function(x) {
   res<-chwirut2.res(x)
   f<-sum(res*res)
}

chwirut2.res <- function(b) {
   xx<-Chwirut2$x # note case !!
   yy<-Chwirut2$y
   res <- rep(NA, length(xx))
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
   res<-exp(-b1*xx)/(b2+b3*xx) - yy
   return(res)
}

# chwirut2 - Jacobian
chwirut2.jac <- function(b) {
xx<-Chwirut2$x # note caps!
J <- matrix(0, length(xx), 3) # define the size of the Jacobian
   b1<-b[1]
   b2<-b[2]
   b3<-b[3]
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

chwirut2.h <- function(x) {
stop("not defined")
   JJ<-chwirut2.jac(x)
   H <- t(JJ) %*% JJ
   res<-chwirut2.res(x)
stop("not defined")

}

chwirut2.g<-function(x) {
#   stop("not defined")
   JJ<-chwirut2.jac(x)
   res<-chwirut2.res(x)
   gg<-as.vector(2.0*t(JJ) %*% res)
   return(gg)
}

chwirut2.fgh<-function(x) {
   f<-chwirut2.f(x)
   g<-chwirut2.g(x)
   H<-chwirut2.h(x)
   fgh<-list(value=f,gradient=g,hessian=H)
}

chwirut2.setup<-function() {
   library(NISTnls) # get parent collection
   data(Chwirut2) # and load up the data into x and y
   start<-c( 0.1, 0.01, 0.02)
}

chwirut2.test<-function() {
}   
