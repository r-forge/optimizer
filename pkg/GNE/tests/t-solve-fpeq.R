require(GNE)
if(FALSE)
{
x0 <- 10 
merit <- function(x, y) 
{
  cat("x", class(x), "\n")
  print(x)
  list(value=sum(exp(x)-1-x)^2+sum(exp(y)-1-y)^2, counts=c(1,0), iter=1)
}
fn <- function(x) 
{
  cat("x", class(x), "\n")
  print(x)
  list(value=exp(x)-1, counts=c(1,0), iter=1)
}

GNE:::fpeq(x0, fn, merit, "RRE", silent=FALSE, control = list(echo=3))
traceback()

mystep <- 3.4
Delta <- 1:3

mystep*Delta
Delta*mystep
Delta*1:2
array(1, 1) + 1:2
as.vector(array(1, 1)) + 1:2
}