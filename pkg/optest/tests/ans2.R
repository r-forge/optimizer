require(graphics)
require(optest)

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
##JN Since this is a single method, the details could have wrong structure
start <- c(-1.2,1)

# methlist <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "SANN", "nlm", "nlminb", 
#               "lbfgsb3", "Rcgmin", "Rtnmin", "Rvmmin", "spg", "ucminf", 
#               "newuoa", "bobyqa", "uobyqa", "nmkb", "hjkb")
# methlist <- c("lbfgsb3")

for (meth in methlist){ 
   myopt <- optest(start, fr, gr=grr, hess=NULL, lower=-Inf, upper=Inf, 
            method=meth, itnmax=NULL, hessian=FALSE, control=list(trace=2))
   cat("myopt for method = ",meth,":\n")
   print(myopt)
   tmp<-readline("continue")
}

