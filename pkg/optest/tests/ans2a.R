# require(graphics)
# require(optest)
rm(list=ls())

fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
    x1 <- x[1]
    x2 <- x[2]
    c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
       200 *      (x2 - x1 * x1))
}

####
require(numDeriv)

# source("optest/R/op.R")
# source("optest/R/optest.setup.R")
source("optest/R/bmchk.R")
source("optest/R/fnchk.R")
source("optest/R/optimr.R")
source("optest/R/scalecheck.R")
source("optest/R/ctrldefault.R")
require(lbfgsb3)
require(Rcgmin)
require(Rvmmin)
require(Rtnmin)
require(minqa)
require(dfoptim)
require(ucminf)
require(BB)

##JN Since this is a single method, the details could have wrong structure
start <- c(-1.2,1)

methlist <- c("lbfgsb3", "Rcgmin", "Rtnmin", "Rvmmin", "spg", "ucminf", 
              "newuoa", "bobyqa", "uobyqa", "nmkb", "hjkb", "BFGS", "CG", "Nelder-Mead", 
               "L-BFGS-B", "SANN", "nlm", "nlminb")

# myenv <- optest.setup(start, fr, grr)
for (meth in methlist){ 
   msg <- paste("Optimr attempt using ",meth)
   tmp<-readline(msg)
#   print(myenv$spar)
#   print(myenv$ufn)
#   print(myenv$ugr)
#   print(meth)
   
   mydo <- optimr(start, fr, grr, method=meth, control=list(trace=1))
   print(mydo)
}
# rm(myenv) # to try to remove duplicate defs


