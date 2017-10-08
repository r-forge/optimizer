## @knitr ##Chwirut1.prb
# This is file ##Chwirut1.prb
probname <- "##Chwirut1"
probdesc <- "Put your description in double quotes.
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
##Chwirut1.formula <- ( y ~ b1*x**b2 )

#- setup

library("NISTnls", character.only=TRUE)
mypdata <- eval(parse(text=data("Chwirut1")))

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

Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Chwirut1)
Try(fm1 <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut1, trace = TRUE,
               start = c(b1 = 0.1, b2 = 0.01, b3 = 0.02)))
Try(fm1a <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut1, trace = TRUE,
                start = c(b1 = 0.1, b2 = 0.01, b3 = 0.02), alg = "port"))
Try(fm2 <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut1, trace = TRUE,
               start = c(b1 = 0.15, b2 = 0.008, b3 = 0.010)))
Try(fm2a <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut1, trace = TRUE,
                start = c(b1 = 0.15, b2 = 0.008, b3 = 0.010), alg = "port"))
Try(fm3 <- nls(y ~ exp(-b1*x)/(1+p3*x), data = Chwirut1, trace = TRUE,
               start = c(b1 = 0.1, p3 = 0.02/0.01), algorithm = "plinear"))
Try(fm4 <- nls(y ~ exp(-b1*x)/(1+p3*x), data = Chwirut1, trace = TRUE,
               start = c(b1 = 0.15, p3 = 0.01/0.008), algorithm = "plinear"))

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
  

