## Watson test
context("Watson test")
test_that("Watson tests", {
    wx0 <- function (n = 15)
    {
        if (!(2 <= n && n <= 31)) {
            stop("Watson: n must be between 2-31")
        }
        rep(0, n)
    }
    waf <- function (par)
    {
        n <- length(par)
        if (!(2 <= n && n <= 31)) {
            stop("Watson: n must be between 2-31")
        }
        fsum <- 0
        for (i in 1:29) {
            ti <- i/29
            fa <- 0
            fb <- 0
            tij2 <- ti^(0:(n - 2))
            tij1 <- ti^(0:(n - 1))
            fa <- sum(1:(n - 1) * par[2:n] * tij2)
            fb <- sum(par * tij1)
            fi <- fa - fb * fb - 1
            fsum <- fsum + fi * fi
        }
        f30 <- par[1]
        f31 <- par[2] - par[1]^2 - 1
        fsum + f30 * f30 + f31 * f31
    }
    wag <- function (par)
    {
        n <- length(par)
        if (!(2 <= n && n <= 31)) {
            stop("Watson: n must be between 2-31")
        }
        grad <- rep(0, n)
        for (i in 1:29) {
            ti <- i/29
            fa <- 0
            fb <- 0
            tij2 <- ti^(-1:(n - 2))
            tij1 <- ti^(0:(n - 1))
            fa <- sum(1:(n - 1) * par[2:n] * tij2[2:length(tij2)])
            fb <- sum(par * tij1)
            fi <- fa - fb * fb - 1
            fi2 <- 2 * fi
            grad <- grad + fi2 * (0:(n - 1) * tij2 - 2 * fb * tij1)
        }
        f31 <- par[2] - par[1]^2 - 1
        grad[1] <- grad[1] + 2 * par[1]
        grad[1] <- grad[1] - 4 * par[1] * f31
        grad[2] <- grad[2] + 2 * f31
        grad
    }
    n <- 12
    x0 <- wx0(n)

    w12lbfsgb3c1 <- lbfgsb3c(x0, waf, wag, control=list(trace=0, maxit=100))
    w12lbfsgb3c <- lbfgsb3c(x0, waf, wag, control=list(trace=0, maxit=200))


})
