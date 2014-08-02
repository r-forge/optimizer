## Optimization test function GENROSE
## ?? refs (put in .doc??)
rm(list = ls())
library(lbfgsb3)

genrose.f <- function(x, gs = NULL) {
    # objective function
    ## One generalization of the Rosenbrock banana valley
    #   function (n parameters)
    n <- length(x)
    if (is.null(gs)) {
        gs = 100
    }
    fval <- 1 + sum(gs * (x[1:(n - 1)]^2 - x[2:n])^2 + (x[2:n] - 
        1)^2)
    return(fval)
}

genrose.g <- function(x, gs = NULL) {
    # vectorized gradient for genrose.f
    # Ravi Varadhan 2009-04-03
    n <- length(x)
    if (is.null(gs)) {
        gs = 100
    }
    gg <- as.vector(rep(0, n))
    tn <- 2:n
    tn1 <- tn - 1
    z1 <- x[tn] - x[tn1]^2
    z2 <- 1 - x[tn]
    gg[tn] <- 2 * (gs * z1 - z2)
    gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
    gg
}



cat("\n\n Unconstrained test\n")
xx <- rep(3, 1000)
lo <- -Inf
up <- Inf
t1000u <- system.time(ans1000u <- lbfgsb3(xx, genrose.f, 
    genrose.g, gs = 10))[1]
cat("final fn value =", ans1000u$value, "\n")
cat("time = ", t1000u, "\n")
t1000uo <- system.time(ao1000u <- optim(xx, genrose.f, 
    genrose.g, method = "L-BFGS-B", gs = 10, control=list(trace=1)))[1]
cat("final fn value =", anso1000u$value, "\n")
cat("time = ", t1000uo, "\n")

t1000un <- system.time(ans1000u <- lbfgsb3(xx, genrose.f, 
    gr = NULL, gs = 10))[1]
cat("final fn value =", ans1000u$value, "\n")
cat("time = ", t1000u, "\n")
t1000uo <- system.time(ao1000u <- optim(xx, genrose.f, 
    gr = NULL, method = "L-BFGS-B", gs = 10, control=list(trace=1)))[1]
cat("final fn value =", anso1000u$value, "\n")
cat("time = ", t1000uo, "\n")
