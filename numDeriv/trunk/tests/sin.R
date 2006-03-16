if(!require("numDeriv"))stop("this test requires numDeriv.")

  # test sin and exp
  func1 <- function(x){ sin(10*x) - exp(-x) }

  curve(func1,from=0,to=5)

  x <- 2.04
  numd1 <- grad(func1, x, method="Richardson")
  exact <- 10*cos(10*x) + exp(-x)
  c(numd1, exact, (numd1 - exact)/exact)

#  numd2 <- grad(func1, x, add options, method="Richardson")
#  exact <- -100*sin(10*x) - exp(-x)
#  c(numd2,exact, (numd1 - exact)/exact)


  x <- c(0:10)
  numd1 <- grad(func1, x, method="Richardson")   # CHECK, NOT WORKING
  exact <- 10*cos(10*x) + exp(-x)
  cbind(numd1, exact, (numd1 - exact)/exact)

#  check jacobian, grad using scalar valued function with vector arg

func2a <- function(x) sum(sin(x))

x <- (0:10)*2*pi/10
#anal.j <- NEED ANALYTIC
#calc.j <- jacobian(func2a, x) # CHECK, NOT WORKING
calc.R <- grad(func2a, x, method="Richardson")
calc.S <- grad(func2a, x, method="simple")

calc.genD <- genD(func2a, x)$D[,1:length(x)]

#if( max(abs(anal.j - calc.R )) > 1e-9) stop("gradient test FAILED")


