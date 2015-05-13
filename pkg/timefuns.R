source("~/R-optimtest/testfnsR/grosenbrock/grose.R")
grose.h <- function(x, gs){
  return(NULL)
}
lfn <- function(x, ...){ # list based function
    fv <- grose.f(x, ...)
    gv <- grose.g(x, ...)
    hv <- grose.h(x, ...)
    fgh <- list(value=fv, gradient=gv, hessia=hv)
}
afn <- function(x, ...){ # attribute based function
    fv <- grose.f(x, ...)
    gv <- grose.g(x, ...)
    hv <- grose.h(x, ...)
    attribute(fv,"gradient") <- gv
    attribute(fv, "hessian") <- hv
    fv
}
require(microbenchmark)
xx<-rep(pi, 10)
microbenchmark(lfn)
microbenchmark(afn)


afgh <- function(x)
{
    gr <- function(x1, x2)
        c(-400*x1*(x2 - x1*x1) - 2*(1-x1), 200*(x2 - x1*x1))
    h <- function(x1, x2) {
        a11 <- 2 - 400*x2 + 1200*x1*x1
        a21 <- -400*x1
        matrix(c(a11, a21, a21, 200), 2, 2)
    }
    x1 <- x[1]; x2 <- x[2]
    res<- 100*(x2 - x1*x1)^2 + (1-x1)^2
    attr(res, "gradient") <- gr(x1, x2)
    attr(res, "hessian") <- h(x1, x2)
    return(res)
}

lfgh <- function(x)
{
    gr <- function(x1, x2)
        c(-400*x1*(x2 - x1*x1) - 2*(1-x1), 200*(x2 - x1*x1))
    h <- function(x1, x2) {
        a11 <- 2 - 400*x2 + 1200*x1*x1
        a21 <- -400*x1
        matrix(c(a11, a21, a21, 200), 2, 2)
    }
    x1 <- x[1]; x2 <- x[2]
    res<- 100*(x2 - x1*x1)^2 + (1-x1)^2
    lfgh <- list(value=res, gradient=gr, hessian=h)
}

microbenchmark(lfgh)
microbenchmark(afgh)


