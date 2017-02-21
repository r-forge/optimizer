# Brown test function
# ?? does it need even n
# ref??

brown.f <- function(x) {
n <- length(x)
## Should we check for odd n and fail out?
if (2*trunc(n/2) != n) { return(NA) }
odd <- seq(1,n,by=2)
even <- seq(2,n,by=2)
sum((x[odd]^2)^(x[even]^2 + 1) + (x[even]^2)^(x[odd]^2 + 1))
}

brown.g <- function(x) {
    stop(" GRADIENT NOT YET DEFINED ")
      
}

brown.setup <- function(n=0) {
    stop(" SETUP NOT YET DEFINED ")
      
}
