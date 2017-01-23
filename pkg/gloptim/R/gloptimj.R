gloptimj <- function(fn, lb, ub, x0 = NULL,
        method = c("deoptim"), type = NULL,
        minimize = TRUE, control = list(), ...) {

    cntrl <- list(info = FALSE,    # shall info/trace be shown
                  popsize = NULL,  # population size
                  itermax = NULL   # max. no. of iterations
                  )

    fun = match.fun(fn)
    f <- function(x) fun(x, ...)

    cat("Global solver/method:", method, "\n")

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
    } else {
        stop("Argument 'method' has not yet been implemented.")
    }
}
