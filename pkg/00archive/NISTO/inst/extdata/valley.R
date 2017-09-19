# valley test function
# Ref: ??


valley.f <- function(x) {
   c1 <- 1.003344481605351
   c2 <- -3.344481605351171e-03
   n <- length(x)
   f <- rep(NA, n)
   j <- 3 * (1:(n/3))
   jm2 <- j - 2
   jm1 <- j - 1
   f[jm2] <- (c2 * x[jm2]^3 + c1 * x[jm2]) * exp(-(x[jm2]^2)/100) - 1
   f[jm1] <- 10 * (sin(x[jm2]) - x[jm1])
   f[j] <- 10 * (cos(x[jm2]) - x[j])
   sum(f*f)
   }

valley.g <- function(x) {
  stop("gradient not yet defined")
}

valley.setup <- function(x) {
  stop("setup not yet defined")

} 
