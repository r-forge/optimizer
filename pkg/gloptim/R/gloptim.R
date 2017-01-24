gloptim <- function(fn, lb, ub, x0 = NULL,
        method = c("deoptim", "deoptimr", "simplede", 
                   "ga",
                   "smco", "soma"),
        type = NULL,
        minimize = TRUE, control = list(), ...) {

    ## Checking input argument for being functions resp. numeric
    fun = match.fun(fn)
    f <- function(x) fun(x, ...)
    
    method = match.arg(method)
    cat("Global solver/method:", method, "\n")

    stopifnot(is.numeric(lb), is.numeric(ub))
    d <- length(lb)
    if (length(ub) != d) stop("Lower and upper bound must be of the same length.")

    if (!is.null(x0)) {
        stopifnot(is.numeric(x0))
        if (length(x0 != d)) stop("")
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
        return(list(xmin = sol$optim$bestmem, fmin = sol$optim$bestval))

    } else if (method == "deoptimr") {
        sol <- DEoptimR::JDEoptim(lower = lb, upper = ub,
                                  fn = fn,
                                  maxiter = cntrl$maxiter)
        return(list(xmin = sol$par, fmin = sol$value))

    } else if (method == "simplede") {
        sol <- adagio::simpleDE(fun = fn, lower = lb, upper = ub,
                                N = cntrl$popsize, nmax = cntrl$maxiter,
                                log = cntrl$info)
        return(list(xmin = sol$xmin, fmin = sol$fmin))

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
    } else if (method == "smco") {
        sol <- smco::smco(par = x0, fn, gr = NULL, ..., N = d, 
                    LB=lb, UB=ub, maxiter = cntrl$maxiter)
        return(list(xmin = sol$par, fmin = sol$value))

    ##
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
