# new broyden tridiagonal to test
broytri2.res <- function(x) {
    n <- length(x)
    tx <- rep(0, n+2)
    # x(i) = tx(i+1)
    # fi = (3-2*x(i))*x(i) - x(i-1)-2*x(i+1)+1
    # x(0)=x(n+1)==0
    xm <- tx[1:n]
    xp <- tx[3:(n+2)]
    f <- (3-2*x)*x -xm -2*xp + 1
}
broytri2.f <- function(x){
  val <- as.numeric(crossprod(broytri2.res(x)))
}
