gloptim <- function(fn, lb, ub, x0 = NULL,
        method = c("deoptim", "ga", "smco", "soma"), type = NULL,
        minimize = TRUE, control = list(), ...) {
    
    fun = match.fun(fn)
    f <- function(x) fun(x, ...)
    
    method = match.arg(method)
    cat("Global solver/method:", method, "\n")

    cntrl <- list(info = FALSE,    # shall info/trace be shown
                  popsize = NULL,  # population size
                  itermax = NULL   # max. no. of iterations
                  )
    for (nm in names(control)) {
        if (nm %in% names(cntrl)) {
            cntrl[nm] <- control[nm]
        } else {
            stop("Unknown name in control list: '", nm, "'.", call. = FALSE)
        }
    }

    if (method == "ga") {
        if (minimize) s <- -1 else s <- 1
        fn <- function(x) s * f(x)
        if (is.null(cntrl$popsize)) popSize <- 100 else popSize <- cntrl$popsize
        if (is.null(cntrl$itermax)) maxiter <- 100 else maxiter <- cntrl$itermax

        sol <- GA::ga(type = "real-valued", fitness = fn,
                      min = lb, max = ub,
                      popSize = popSize,
                      maxiter = maxiter,
                      monitor = cntrl$info)
        return(list(xmin = sol@solution,
                    fmin = s * sol@fitnessValue))

    } else if (method == "deoptim") {
        if (is.null(cntrl$itermax)) maxiter <- 1000 else maxiter <- cntrl$itermax
        
        sol <- DEoptim::DEoptim(fn, lower = lb, upper = ub,
                                DEoptim::DEoptim.control(
                                trace = cntrl$info, itermax = cntrl$itermax))
        return(list(xmin = sol$optim$bestmem, fmin = sol$optim$bestval))
  } else if (method == "smco") {
    if (is.null(cntrl$itermax)) maxiter <- 1000 else maxiter <- cntrl$itermax
    if (is.numeric(x0)) N <- length(x0) else N <- max(length(lb), length(ub))
    ## ?? improve this
    sol <- smco(par = x0, fn, gr = NULL, ..., N = N, 
                LB=lb, UB=ub, maxiter = maxiter)
    
    return(list(xmin = sol$optim$bestmem, fmin = sol$optim$bestval))    
  } else if (method == "soma") {
    if (is.null(cntrl$itermax)) maxiter <- 1000 else maxiter <- cntrl$itermax
    sol <- soma(costFunction = fn, bounds=list(min=lb, max=ub))
                
    return(list(xmin = sol$optim$bestmem, fmin = sol$optim$bestval))    
  } else {
        stop("Argument 'method' has not yet been implemented.")
   }
}
