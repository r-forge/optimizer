## @knitr Nesterov.prb
# This is file Nesterov.prb
probname <- "Nesterov"
probdesc <- "Put your description in double quotes.
"

#- Note: environment / list "pe" must already exist

if (! exists("pe")) {stop("Environment 'pe' must exist before running problems.")}

#- nls format expression
Nesterov.formula <- ( y ~ b1*x**b2 )

#- setup
# Nesterov.setup<-function() {
y <- c( 2.138, 3.421, 3.597, 4.340, 4.882, 5.660) # example only
x <- c( 1.309, 1.471, 1.490, 1.565, 1.611, 1.680)
Nesterov.df <- data.frame(x,y)
rm(x)
rm(y)

Nesterov.start <- function(indx) {
  ...

  start
}
#- Problem types will be determined from available functions / formulas
#- ?? We could have problems with quotation marks WITHIN the code.
#- lower = NA # if bounds are present, then we use bounded methods
#- upper = NA # 

#- function
Nesterov.f <- function(x) { # ?? should remove -- this is a sumsquares problem
    n <- length(x)
    f <- (1 - x[1])^2/4
    for (i in 1:(n - 1)) {
        f <- f + (1 + x[i + 1] - 2 * x[i]^2)^2
    }
    f
}

#- gradient
Nesterov.g<-function(x) {
    n <- length(x)
    g <- rep(0, n)
    g[1] <- (x[1] - 1)/2
    for (i in 1:(n - 1)) {
        r = 1 + x[i + 1] - 2 * x[i]^2
        g[i + 1] <- g[i + 1] + 2 * r
        g[i] <- g[i] - 8 * x[i] * r
    }
    g
}

#- hessian
Nesterov.h <- function(x) {
#- THIS IS NOT COMPLETE??
  stop("DanielWood.h IS NOT COMPLETE")
  res<-DanielWood.res(x)
  JJ<-DanielWood.jac(x)
  pe$khess <- pe$khess + 1
  H <- t(JJ) %*% JJ
}


#- residual
Nesterov.res <- function(b) {
return(res)
}

#- Jacobian
# DanielWood - Jacobian
Nesterov.jac <- function(b) {
return(J)
}

#- nleq -- ?? section if appropriate


#- test
## put example calls of the function, possibly including calls to 
# optimizations and nonlinear least squares etc.



#- End DanielWood.prb   

