##################################################################
cgdefault <- function(npar) { 

# offset changed from 100 to 1000 on 180710
      ctrl.default <- list(
        acctol = 0.0001, # used for acceptable point test in backtrack linesearch
      	badval = (0.5)*.Machine$double.xmax, # use this value as a flag of non-computable item
        bigval = .Machine$double.xmax*0.01, # a very large number (note smaller than badval)
        defgrapprox = "grfwd", # use forward approximation as default. Could argue for grcentral
        dowarn = TRUE,  # generally want to turn on warnings
        eps = 1e-07, # a tolerance for a small quantity (about single precision level)
        epstol = .Machine$double.eps, # but this is the machine epsilon
        have.bounds = FALSE, # normally have UNCONSTRAINED function
        maximize = FALSE, # normally MINIMIZE (see fnscale)
        maxit = 500*round(sqrt(npar+1)), # limit on number of iterations or gradient evaluations
        maxfeval = 500*round(sqrt(npar+1)), # limit on function evaluations
        newtstep=1, # default stepsize is 1 (Newton stepsize)
        offset = 1000.0, # used for equality test (a + offset) == (b + offset)
        stepinc = 10,
        steplen0 = 0.75,
        steplenn = 1.0, # ?? fix choices of names 
        stepredn = 0.2,
        tol =  npar * (npar * .Machine$double.eps),  # for gradient test.
        # Note -- integer overflow if npar*npar*.Machine$double.eps, 
        trace = 0,
        cgoffset = 100,
        cgstepredn = 0.15,
        cgoldstep = 0.8,
        cgstinflate = 2.25,
        cgstep0max = 1,
        cgfminplus = 0 # how many abs(fmin) to add to fmin allowed for quad inv interp??
      )
}
