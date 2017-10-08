## @knitr ##Chwirut2.prb
# This is file ##Chwirut2.prb
probname <- "##Chwirut2"
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
##Chwirut2.formula <- ( y ~ b1*x**b2 )

#- setup

library("NISTnls", character.only=TRUE)
mypdata <- eval(parse(text=data("Chwirut2")))

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

## Examples

Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = Chwirut2)
Try(fm1 <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut2, trace = TRUE,
               start = c(b1 = 0.1 , b2 = 0.01, b3 = 0.02)))
Try(fm1a <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut2, trace = TRUE,
                start = c(b1 = 0.1 , b2 = 0.01, b3 = 0.02), alg = "port"))
Try(fm2 <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut2, trace = TRUE,
               start = c(b1 = 0.15 , b2 = 0.008, b3 = 0.01)))
Try(fm2a <- nls(y ~ exp(-b1*x)/(b2+b3*x), data = Chwirut2, trace = TRUE,
                start = c(b1 = 0.15 , b2 = 0.008, b3 = 0.01), alg = "port"))
Try(fm3 <- nls(y ~ exp(-b1*x)/(1+p3*x), data = Chwirut2, trace = TRUE,
               start = c(b1 = 0.1, p3 = 2.), alg = "plinear"))
Try(fm4 <- nls(y ~ exp(-b1*x)/(1+p3*x), data = Chwirut2, trace = TRUE,
               start = c(b1 = 0.15, p3 = 0.01/0.008), alg = "plinear"))
