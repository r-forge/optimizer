fnHald <- function(x) {
    stopifnot(is.numeric(x), length(x) == 5)
    if (all(x == 1)) return(exp(1))
    t <- -1 + (c(1:21) - 1)/10
    v <- (x[1]+x[2]*t) / (1 + x[3]*t + x[4]*t^2 + x[5]*t^3) - exp(t)
    max(abs(v))
  }  

library(optimrx)

st0 <- seq(1:5)/10

sol0unc <- opm(st0, fnHald, gr="grfwd", method="ALL", control=list(trace=1))
summary(sol0unc, order=value)

sol0uncnogr <- opm(st0, fnHald, gr=NULL, method="ALL", control=list(trace=1))
summary(sol0uncnogr, order=value)


