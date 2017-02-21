fnNest<-function(x){
    n <- length(x)
    f <- (1 - x[1])^2/4
    for (i in 1:(n - 1)) {
        f <- f + (1 + x[i + 1] - 2 * x[i]^2)^2
    }
    f
}

grNest <- function(x) {
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
