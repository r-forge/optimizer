##
##  g l o p t i m . R
##


gloptim <- function(fn, lb, ub,
                    x0 = NULL, rand = FALSE,
        method = c("deoptim", "cppdeoptim", "deoptimr",         # **DE**
                   "deopt", "simplede", "simpleea",             # **EA**
                   "gensa", "ga", "genoud", # "rbga"            # **GA**
                   "pso", "psopt", "hydropso", # "psoptim"      # **PSO**
                   "direct", "crs2lm", "isres", "stogo",        # **NLoptr**
                   "cmaoptim", "cmaes", "cmaesr", "purecmaes",  # **CMA-ES**
                   "malschains", "ceimopt",                     # **CE**
                   "smco", "soma"), # "tabusearch"              # --others--
        type = NULL,
        g = NULL, gr = NULL,
        minimize = TRUE, control = list(), ...) {

    ## Checking input argument for being functions resp. numeric
    fun = match.fun(fn)
    f <- function(x) fun(x, ...)

    method = match.arg(method)
    cat("Global solver/method:", method, "\n")

    stopifnot(is.numeric(lb), is.numeric(ub))
    d <- length(lb)
    if (length(ub) != d)
        stop("Lower and upper bound must be of the same length.")

    if (!is.null(x0)) {
        stopifnot(is.numeric(x0))
        if (length(x0) != d)
            stop("Argument 'x0' must have the same length as bounds.")
        if (any(x0 < lb) || any(x0 > ub))
            stop("Argument 'x0' not between lower and upper bound.")
    } else {
        if (rand) {
            x0 <- lb + runif(d) * (ub - lb)
        } else {
            x0 <- (ub + lb) / 2.0
        }
    }

    ## Handling of requested method and control options 
    cntrl <- list(
                  popsize =  20*d,      # population size
                  maxiter = 200*d,      # max. no. of iterations
                  info    = FALSE       # shall info/trace be shown
    )
    for (nm in names(control)) {
        if (nm %in% names(cntrl)) {
            cntrl[nm] <- control[nm]
		} else
		    stop("Unknown name in control list: '", nm, "'.", call. = FALSE)
    }

    ## Differential Evolution (DE) in packages DEoptim and DEoptimR
    if (method == "deoptim") {
        sol <- DEoptim::DEoptim(fn, lower = lb, upper = ub,
                                DEoptim::DEoptim.control(
                                trace = cntrl$info, itermax = cntrl$maxiter))
        return(list(xmin = sol$optim$bestmem, fmin = sol$optim$bestval,
                    niter = sol$optim$iter, nfeval = sol$optim$nfeval,
                    commment = ""))

    } else if (method == "cppdeoptim") {
        sol <- RcppDE::DEoptim(fn, lower = lb, upper = ub,
                                RcppDE::DEoptim.control(
                                trace = cntrl$info, itermax = cntrl$maxiter))
        return(list(xmin = sol$optim$bestmem, fmin = sol$optim$bestval,
                    niter = sol$optim$iter, nfeval = sol$optim$nfeval,
                    commment = ""))
        

    } else if (method == "deoptimr") {
        sol <- DEoptimR::JDEoptim(lower = lb, upper = ub,
                                  fn = fn,
                                  maxiter = cntrl$maxiter)
        return(list(xmin = sol$par, fmin = sol$value,
                    niter = sol$iter, nfeval = NA,
                    comment = ""))

    } else if (method == "deopt") {
        sol <- NMOF::DEopt(OF = fn,
                           algo = list(nP = cntrl$popsize, nG = cntrl$maxiter,
                                       min = lb, max = ub,
                                       printDetail = 2 * cntrl$info,
                                       printBar = FALSE))
        return(list(xmin = sol$xbest, fmin = sol$OFvalue,
                    niter = cntrl$maxiter, nfeval = NA,
                    comment = ""))

    } else if (method == "simplede") {
        sol <- adagio::simpleDE(fun = fn, lower = lb, upper = ub,
                                N = cntrl$popsize, nmax = cntrl$maxiter,
                                log = cntrl$info)
        return(list(xmin = sol$xmin, fmin = sol$fmin))

    } else if (method == "simpleea") {
        sol <- adagio::simpleEA(fn = fn, lower = lb, upper = ub,
                                # N = cntrl$popsize,
                                log = cntrl$info)
        return(list(xmin = sol$par, fmin = sol$val))
        
    ## Simulated Annealing (SA) in package GenSA
    } else if (method == "gensa") {
        sol <- GenSA::GenSA(fn=fn, lower=lb, upper=ub,
                            control=list(maxit=cntrl$maxiter,
                                         smooth=FALSE,
                                         trace.mat=cntrl$info))
        return(list(xmin = sol$par, fmin = sol$value,
                    niter = cntrl$maxiter, feval = sol$counts))

    ## Genetic Algorithms (GA) in package GA
    } else if (method == "ga") {
        if (minimize) s <- -1 else s <- 1   # Default in GA is maximization
        fn <- function(x) s * f(x)
        sol <- GA::ga(type = "real-valued", fitness = fn,
                      min = lb, max = ub,
                      popSize = cntrl$popsize,
                      maxiter = cntrl$maxiter,
                      monitor = cntrl$info)
        return(list(xmin = sol@solution, fmin = s * sol@fitnessValue))

    ## 
    } else if (method == "genoud") {
        if (is.null(gr))
            gr <- function(x) rep(0, d)
        sol <- rgenoud::genoud(fn = fn, nvars = d, gr = gr,
                               # pop.size = cntrl$popsize,
                               # max.generations = cntrl$maxiter,
                               Domains = cbind(lb, ub), boundary.enforcement = 1,
                               print.level = 0)
        return(list(xmin = sol$par, fmin = sol$value))

    ## Particle Swarm Optimization (PSO) in packages pso, psoptim, and NMOF
    } else if (method == "pso") {
        if (is.null(x0)) x0 <- rep(NA, d)
        sol <- pso::psoptim(par=x0, fn = fn,
                            lower = lb, upper = ub,
                            control=list(maxit=cntrl$maxiter))  # maxf=14*maxit
        return(list(xmin = sol$par, fmin = sol$value))

    } else if (method == "psopt") {
        sol <- NMOF::PSopt(OF = fn,
                           algo = list(nP = cntrl$popsize, nG = cntrl$maxiter,
                                       min = lb, max = ub,
                                       printDetail = 0, printBar = FALSE))
        return(list(xmin = sol$xbest, fmin = sol$OFvalue))

## Package psoptim, recurring error:
## Error in plot.window(...) : need finite 'ylim' values
#     } else if (method == "psoptim") {
#         fnew <- function(x) {
#             if (is.vector(x)) {
#                 fval <- fn(x)
#             } else if (is.matrix(x)) {
#                 fval <- numeric(nrow(x))
#                 for (i in 1:nrow(x)) fval[i] <- fn(x[i, ])
#             }
#             return(fval)
#         }
#         sol <- psoptim::psoptim(FUN = fnew,
#                                 n = cntrl$popsize, max.loop = cntrl$maxiter,
#                                 xmin = lb, xmax = ub, vmax = rep(1, d),
#                                 anim = FALSE)
#         return(list(xmin = sol$sol, fmin = sol$val))

    } else if (method == "hydropso") {
        sol <- hydroPSO::hydroPSO(par = x0, fn = fn,
                                  lower = lb, upper = ub,
                                  control = list(npart = cntrl$popsize,
                                                 maxit = cntrl$maxiter,
                                                 verbose = cntrl$info))
        return(list(xmin = sol$par, fmin = sol$value))

    ## Cross Entropy (CE) inspired methods in packages RCEIM and CEoptim
    } else if (method == "ceimopt") {
        sol <- RCEIM::ceimOpt(OptimFunction = fn, nParams = d,
                              Ntot = cntrl$popsize, maxIter = cntrl$maxiter,
                              boundaries = cbind(lb, ub),
                              verbose = cntrl$info)
        xmin <- unname(sol$BestMember[1:d])
        fmin <- unname(sol$BestMember[d+1])
        return(list(xmin = xmin, fmin = fmin))

    ## Different stochastic solvers in package nloptr
    } else if (method == "direct") {
        sol <- nloptr::direct(fn = fn, lower = lb, upper = ub,
                              control = list(maxeval = 10*cntrl$maxiter),
                              nl.info = cntrl$info)
        return(list(xmin = sol$par, fmin = sol$value))

    } else if (method == "crs2lm") {
        sol <- nloptr::crs2lm(x0 = x0, fn = fn,
                              lower = lb, upper = ub,
                              maxeval = 10*cntrl$maxiter, pop.size = cntrl$popsize,
                              nl.info = cntrl$info)
        return(list(xmin = sol$par, fmin = sol$value))

    } else if (method == "isres") {
        sol <- nloptr::isres(x0 = x0, fn = fn,
                             lower = lb, upper = ub,
                             maxeval = 10*cntrl$maxiter, pop.size = cntrl$popsize,
                             xtol_rel = 1e-7, nl.info = cntrl$info)
        return(list(xmin = sol$par, fmin = sol$value))

    } else if (method == "stogo") {
        sol <- nloptr::stogo(x0 = x0, fn = fn, gr = NULL,
                             lower = lb, upper = ub,
                             maxeval = cntrl$maxiter,
                             nl.info = cntrl$info)
        return(list(xmin = sol$par, fmin = sol$value))

    ## Covariance Matrix Adaptation Evolution Strategy(CMA-ES)
    ## in packages rCMA, parma, cmaes, adagio, and Rmalschains
    } else if (method == "cmaoptim") {
        cma_obj <- rCMA::cmaNew()
        rCMA::cmaInit(cma_obj, dimension=d)
        sol <- rCMA::cmaOptimDP(cma_obj, fitFunc = fn,
                                verbose = 0)
        return(list(xmin = sol$bestX, fmin = sol$bestFitness))

    } else if (method == "cmaes") {
        cntrl_cmaes <- parma::cmaes.control()
        cntrl_cmaes$options$MaxIter <- cntrl$maxiter
        cntrl_cmaes$options$PopSize <- cntrl$popsize
        cntrl_cmaes$options$DispModulo <- 0
        cntrl_cmaes$options$DispFinal <- FALSE
        sol <- parma::cmaes(pars = x0, fun = fn,
                            lower = lb, upper = ub,
                            ctrl = cntrl_cmaes)
        return(list(xmin = sol$bestever$x, fmin = sol$bestever$f))

    } else if (method == "cmaesr") {
        fn_parset <- 
            ParamHelpers::makeNumericParamSet("R", len=d,
                                lower=lb, upper=ub)
        fn_fun <- 
            smoof::makeSingleObjectiveFunction("Function", fn = fn,
                                        par.set = fn_parset)
        sol <- cmaesr::cmaes(fn_fun, start.point = x0,
                             monitor = NULL)
        return(list(xmin = sol$best.param, fmin = sol$best.fitness))

    } else if (method == "purecmaes") {
        sol <- adagio::pureCMAES(par = x0, fun = fn,
                                 lower = lb, upper = ub,
                                 stopeval=5*d*cntrl$maxiter)
        return(list(xmin = sol$xmin, fmin = sol$fmin))

    } else if (method == "malschains") {
        sol <- Rmalschains::malschains(fn = fn, dim = d,
                                       lower = lb, upper = ub,
                                       maxEvals = cntrl$popsize * cntrl$maxiter,
                                       verbosity = 0)
        return(list(xmin = sol$sol, fmin = sol$fitness))

    ## Monte Carlo Optimizer (MCO) in package smco
    } else if (method == "smco") {
        sol <- smco::smco(fn = fn, gr = NULL, N = d, 
                    LB = lb, UB = ub,
                    maxiter = cntrl$maxiter, trc = cntrl$info)
        return(list(xmin = sol$par, fmin = sol$value))

    ## Self-Organizing Optimization (SOM) algorithm in package soma
    } else if (method == "soma") {
        sol <- soma::soma(costFunction = fn,
                    bounds=list(min=lb, max=ub) )
        best <- sol$leader
        xmin <- (sol$population)[ ,best]
        return(list(xmin = xmin, fmin = (sol$cost)[best]) )

    } else {
        stop("Argument '", method, "' has not (yet) been implemented.")
    }

} # EoF
