try## rosen2 -- standard 2 parameter Rosenbrock banana shaped valley

require(optimx)

XRosenbrock.f <- function (x) 
{
    n <- length(x)
    x1 <- x[2:n]
    x2 <- x[1:(n - 1)]
    sum(100 * (x1 - x2^2)^2 + (1 - x2)^2)
}
XRosenbrock.g <- function (x) 
{
    n <- length(x)
    g <- rep(NA, n)
    g[1] <- 2 * (x[1] - 1) + 400 * x[1] * (x[1]^2 - x[2])
    if (n > 2) {
        ii <- 2:(n - 1)
        g[ii] <- 2 * (x[ii] - 1) + 400 * x[ii] * (x[ii]^2 - x[ii + 
            1]) + 200 * (x[ii] - x[ii - 1]^2)
    }
    g[n] <- 200 * (x[n] - x[n - 1]^2)
    g
}

xstart <- c(-1.2,1)


ansone<-optimr(xstart, XRosenbrock.f, XRosenbrock.g, method="Rvmmin")

ansone.target <- structure(list(par = c(1, 1), value = 0, counts = structure(c(59, 
39), .Names = c("function", "gradient")), convergence = 2, 
message = "Rvmminu appears to have converged"), .Names = c("par", 
"value", "counts", "convergence", "message"))

ansall <- opm(xstart, XRosenbrock.f, XRosenbrock.g, method="ALL", control=list(kkt=FALSE))
ansall.sum <- summary(ansall, order=value)
print(ansall.sum)

ansallx <- optimx(xstart, XRosenbrock.f, XRosenbrock.g, control=list(all.methods=TRUE, kkt=FALSE))
ansallx.sum <- summary(ansallx, order=value)
print(ansallx.sum)

# test for alternate call of all methods
# ansalla <- opm(xstart, XRosenbrock.f, XRosenbrock.g, control=list(all.methods=TRUE, kkt=FALSE))
# ansalla.sum <- summary(ansalla, order=value)
# print(ansalla.sum)
