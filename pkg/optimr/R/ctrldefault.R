##################################################################
ctrldefault <- function(npar) { 
#  THIS VERSION FOR optimr (reduced set of methods)
#
     ## These are DEFAULTS. They may be nonsense in some contexts.

      allmeth <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "nlm", "nlminb", 
               "Rcgmin", "Rvmmin", "hjn")



      bdmeth <- c("L-BFGS-B", "nlminb", "Rcgmin", "Rvmmin",  "hjn")


      maskmeth <- c("Rcgmin", "Rvmmin", "hjn")
 
      ctrl.default <- list(
        acctol = 0.0001, 
        all.methods = FALSE,
        allmeth = allmeth,
	badval = (0.5)*.Machine$double.xmax,
        bdmeth = bdmeth,
        defgrapprox = "grfwd",
        dowarn = TRUE, 
        eps = 1e-07, 
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
        parchanged = FALSE, 
        parscale = rep(1, npar),
        reltest = 100.0,
        save.failures = TRUE,
	scaletol = 3, 
        steplen0 = 0.75, 
        stepredn = 0.2,
        stopbadupdate = FALSE,
        tol = 0, 
	trace = 0
      )
}
##################################################################

