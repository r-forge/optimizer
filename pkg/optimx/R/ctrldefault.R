##################################################################
ctrldefault <- function(npar) { 
# THIS IS FULL VERSION FOR optimrx
#
     ## These are DEFAULTS. They may be nonsense in some contexts.

      allmeth <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "nlm", "nlminb", 
                "lbfgsb3", "Rcgmin", "Rtnmin", "Rvmmin", "snewton", "snewtonm",
                 "spg", "ucminf", "newuoa", "bobyqa", "nmkb", "hjkb", "hjn", 
                 "lbfgs", "subplex")

#  allpkg has package where element of allmeth is found
      allpkg <-  c("stats", "stats", "stats", "stats", "stats", "stats",
                "lbfgsb3", "Rcgmin", "Rtnmin", "Rvmmin", "optimx", "optimx",
                "BB", "ucminf", "minqa", "minqa", "dfoptim", "dfoptim", 
                "optimx", "lbfgs", "subplex")

     # 160628: uobyqa removed as it fails hobbs from 1,1,1 unscaled

      bdmeth <- c("L-BFGS-B", "nlminb", "lbfgsb3", "Rcgmin", "Rtnmin", "Rvmmin",  
                "bobyqa", "nmkb", "hjkb", "hjn")

      maskmeth <- c("Rcgmin", "Rvmmin", "hjn")
 
      ctrl.default <- list(
        acctol = 0.0001, 
        all.methods = FALSE,
        allmeth = allmeth,
        allpkg = allpkg,
      	badval = (0.5)*.Machine$double.xmax,
        bdmeth = bdmeth,
        bigval = .Machine$double.xmax*0.01,
        defgrapprox = "grfwd",
        defmethod = "Nelder-Mead",
        defstep=1,
        dowarn = TRUE, 
        eps = 1e-07, 
        epstol = .Machine$double.eps,
        fnscale = 1.0, 
        grtesttol=(.Machine$double.eps)^(1/3), 
        have.bounds = FALSE,
        hessasymtol = 0.0001,
        hesstesttol=(.Machine$double.eps)^(1/3), 
        keepinputpar = FALSE,
	      kkt = TRUE,
	      kkttol = 0.001,
	      kkt2tol = 1.0E-6,
        maskmeth = maskmeth,
	      maximize = FALSE,
        maxit = 500*round(sqrt(npar+1)),
	      maxfeval = 5000*round(sqrt(npar+1)),
        offset = 100.0,
        parchanged = FALSE, 
        parscale = rep(1, npar),
        reltest = 100.0,
        save.failures = TRUE,
      	scaletol = 3, 
        stepdec = 0.2, 
        steplen0 = 0.75, 
        stepmax = 5,
        stepmin = 0,
        stepredn = 0.2,
        stopbadupdate = FALSE,
        tol = 0, 
	      trace = 0,
	      watch = FALSE
      )
}
##################################################################

