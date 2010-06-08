##' Handle the common arguments to all the minimizers
##'
##' Establish defaults for the elements of the control list and parse
##' the list that was provided.  The list is converted to an
##' environment.
##'
##' @param ctrl list of control settings
##' @param n length of the par vector
##' 
##' @return an environment containing the control settings
commonArgs <- function(par, fn, ctrl, rho) {
    rho$n <- n <- length(rho$par <- as.double(par))
    stopifnot(all(is.finite(par)),# n > 1,
              is.function(fn),
              length(formals(fn)) >= 1)
    rho$.par. <- numeric(n)             # argument for internal function
    rho$.feval. <- integer(1)           # function evaluation counter

    ## We use all possible control settings in the default.
    ## Extra control settings are ignored.
    cc <- do.call(function(npt = min(n+6L, 2L * n + 1L), rhobeg = NA,
                           rhoend = NA, iprint = 0L, maxfun=10000L,
                           obstop=TRUE, force.start=FALSE)
                  list(npt = npt, rhobeg = rhobeg, rhoend = rhoend,
                       iprint = iprint, maxfun = maxfun, obstop = obstop,
                       force.start = force.start), ctrl)

    ## Create and populate and environment
    ctrl <- new.env(parent = emptyenv()) # ctrl environment should not chain
    lapply(names(cc), function(nm) assign(nm, cc[[nm]], envir = ctrl))

    ## Adjust and check npt
    ctrl$npt <- as.integer(max(n + 2, min(ctrl$npt, (n+1)*(n+2)/2)))
    if (ctrl$npt > (2 * n + 1))
        warning("Setting 'npt' larger than 2 * length(par)+1 not recommended.")

    ## Check and adjust rhobeg and rhoend
    if (is.na(ctrl$rhobeg))
        ctrl$rhobeg <- min(0.95, 0.2 * max(abs(par)))
    if (is.na(ctrl$rhoend)) ctrl$rhoend <- 1.0e-6 * ctrl$rhobeg
    stopifnot(0 < ctrl$rhoend, ctrl$rhoend <= ctrl$rhobeg)

    ## Check recommended range of maxfun
    if (ctrl$maxfun < 10 * n^2)
        warning("'maxfun' less than 10*length(par)^2 not recommended.")
    ctrl
}

##' Nonlinear optimization with box constraints
##'
##' Minimize a function of many variables subject to box constraints
##' by a trust region method ##' that forms quadratic models by 
##' interpolation, using the BOBYQA code written by Mike Powell.
##' 
##' @param par numeric vector of starting parameters (length > 1)
##' @param fn function to be minimized.  The first argument must be
##'     the parameters.
##' @param lower a numeric vector of lower bounds.  If of length 1 it
##'     is expanded.
##' @param upper a numeric vector of upper bounds.  Also may be scalar.
##' @param control a list of control settings
##' @param ... optional, additional arguments to fn
##'
##' @return a list with S3 class bobyqa
##' 
bobyqa <- function(par, fn, lower = -Inf, upper = Inf, control = list(), ...)
{
    ## the "+ 0" is to force a copy of par in the environment.
    ctrl <- commonArgs(par + 0, fn, control, environment())
    n <- length(par)
    fn1 <- function(x) fn(x, ...) # fn1 takes exactly 1 argument
    
    ## check the upper and lower arguments, adjusting if necessary
    lower <- as.double(lower); upper <- as.double(upper)
    if (length(lower) == 1) lower <- rep(lower, n)
    if (length(upper) == 1) upper <- rep(upper, n)
    stopifnot(length(lower) == n, length(upper) == n, all(lower < upper))
    if (any(par < lower | par > upper)) {
        if (ctrl$obstop)
            stop("Starting values violate bounds")
        else {
            par <- pmax(lower, pmax(par, upper))
            warning("Some parameters adjusted to nearest bound")
        }
    }
    rng <- upper - lower

    verb <- 1 < (ctrl$iprint <- as.integer(ctrl$iprint))
    if (verb) {
        cat("npt =", ctrl$npt, ", n = ",n,"\n")
        cat("rhobeg = ", ctrl$rhobeg,", rhoend = ", ctrl$rhoend, "\n")
    }
    if (any(rng < 2 * ctrl$rhobeg)) {
        warning("All upper - lower must be >= 2*rhobeg. Changing rhobeg") 
        rhobeg <- 0.2 * min(rng)
    }
    
    ## Modifications to par if too close to boundary
    if (all(is.finite(upper)) && all(is.finite(lower)) &&
        all(par >= lower) && all(par <= upper) ) {
        if (verb) cat("ctrl$force.start = ", ctrl$force.start,"\n")
        if (!ctrl$force.start) {
            i <- rng < ctrl$rhobeg # Jens modification
            if (any(i)) {
                par[i] <- lower[i] + ctrl$rhobeg
                warning("Some parameters adjusted away from lower bound")
            }
            i <- rng < ctrl$rhobeg      # Jens modification
            if (any(i)) {
                par[i] <- upper[i] - ctrl$rhobeg
                warning("Some parameters adjusted away from upper bound")
            }
        }
    }

    ## force one evaluation of fn1 to check that it works
    fn1(.par.)

    .Call(bobyqa_cpp, par, lower, upper, ctrl, fn1)
}

##' An R interface to the NEWUOA implementation of Powell
##'
##' Minimize a function of many variables by a trust region method
##' that forms quadratic models by interpolation, using the NEWUOA
##' code written by Mike Powell.
##' 
##' @param par numeric vector of starting parameters (length > 1)
##' @param fn function to be minimized.  The first argument must be
##'     the parameters.
##' @param control a list of control settings
##' @param ... optional, additional arguments to fn
##'
##' @return a list with S3 class c("newuoa", "minqa")
newuoa <- function(par, fn, control = list(), ...)
{
    ctrl <- commonArgs(par + 0, fn, control, environment())
    fn1 <- function(x) fn(x, ...)
    fn1(.par.)

    .Call(newuoa_cpp, par, ctrl, fn1)
}

##' An R interface to the UOBYQA implementation of Powell
##'
##' Minimize a function of many variables by a trust region method
##' that forms quadratic models by interpolation, using the UOBYQA
##' code written by Mike Powell.
##' 
##' @param par numeric vector of starting parameters (length > 1)
##' @param fn function to be minimized.  The first argument must be
##'     the parameters.
##' @param control a list of control settings
##' @param ... optional, additional arguments to fn
##'
##' @return a list with S3 class uobyqa
uobyqa <- function(par, fn, control = list(), ...)
{
    ctrl <- commonArgs(par + 0, fn, control, environment())
    fn1 <- function(x) fn(x, ...)
    fn1(.par.)

    .Call(uobyqa_cpp, par, ctrl, fn1)
}

##' Print method for minqa objects (S3)
##'
##' @param x an object of class that inherits from minqa
##' @param digits number of significant digits - doesn't seem to be used
##' @param ... optional arguments.  None are used.
##' @return invisible(x) - side effect is to print
##' @author Douglas Bates
print.minqa <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("parameter estimates:", toString(x$par), "\n")
  cat("objective:", toString(x$fval), "\n")
  cat("number of function evaluations:", toString(x$feval), "\n")
  
  invisible(x)
}
