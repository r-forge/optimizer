# require(graphics)
# require(optest)

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

source("optest/R/op.R")
source("optest/R/optest.setup.R")
source("optest/R/bmchk.R")
source("optest/R/fnchk.R")
source("optest/R/optimr.R")
source("optest/R/scalecheck.R")
source("optest/R/ctrldefault.R")

##JN Since this is a single method, the details could have wrong structure
start <- c(-1.2,1)

methlist <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "SANN", "nlm", "nlminb", 
              "lbfgsb3", "Rcgmin", "Rtnmin", "Rvmmin", "spg", "ucminf", 
              "newuoa", "bobyqa", "uobyqa", "nmkb", "hjkb")

mycfg <- optest.setup(start, fr, grr)

for (meth in methlist){ 
   myopt <- op(start, fr, gr=grr, hess=NULL, 
            method=meth, itnmax=NULL, hessian=FALSE, control=list(trace=1, starttests=FALSE))
   cat("myopt for method = ",meth,":\n")
   print(myopt)
   tmp<-readline("continue to call using mycfg")
   mydo <- optimr(mycfg$spar, mycfg$ufn, mycfg$ugr, method=meth, control=list(trace=2),
        pscale=rep(1,2))
   print(mydo)
   
}

