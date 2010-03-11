bobyqa <- function(par, fn, lower=-Inf, upper=Inf,
                   control = list(), ...)
{
    n <- length(par <- as.double(par))
    if (n < 2)
        stop("bobyqa is not for optimization of single parameter.")

    ## initialize the control list
    ctrl <- list(npt = min(n+6, 2*n+1), rhobeg = NA, rhoend = NA,
                 iprint = 0L, maxfun=10000L, obstop=TRUE,
                 force.start=FALSE)
    ## incorporate any control arguments using partial matches on names
    if (!missing(control)) {
        control <- as.list(control)
        mnms <- pmatch(names(control), names(ctrl)) # matched names
        if (any(is.na(mnms))) {
            warning("Unmatched names:", names(control)[is.na(mnms)],
                    " dropped from control list")
            control <- control[!is.na(mnms)]
            mnms <- mnms[!is.na(mnms)]
            names(control) <- names(ctrl)[mnms]
        }
        ctrl[names(control)] <- control
    }

    ## check the upper and lower arguments, adjusting if necessary
    lower <- as.double(lower); upper <- as.double(upper)
    if (length(lower) == 1) lower <- rep(lower, n)
    if (length(upper) == 1) upper <- rep(upper, n)
    if (length(lower) != n || length(upper) != n)
        stop("Bounds are the wrong length.")
    if (any(upper <= lower)) stop("Overlapping bounds")

    if (any(par < lower | par > upper)) {
        if (ctrl$obstop)
            stop("Starting values violate bounds")
        else {
            par <- pmax(lower, pmax(par, upper))
            warning("Some parameters adjusted to nearest bound")
        }
    }
    rng <- upper - lower
    
    ## Adjust and check npt
    ctrl$npt <- max(n+2, min(as.integer(ctrl$npt), (n+1)*(n+2)/2))
    if (ctrl$npt > (2*n + 1))
        warning("Setting 'npt' larger than 2 * length(par)+1 not recommended.")
    verb <- 1 < (ctrl$iprint <- as.integer(ctrl$iprint))
    if (verb)
        cat("npt is set at ", ctrl$npt, " and n = ",n,"\n")

    ## Check and adjust rhobeg and rhoend
    if (is.na(ctrl$rhobeg))
        ctrl$rhobeg <-
            min(0.95,
                0.2 * ifelse(any(is.finite(rng)), min(rng), max(abs(par))))
    if (is.na(ctrl$rhoend)) {
        ctrl$rhoend <- 1.0e-6 * ctrl$rhobeg
    }
    stopifnot(0 < ctrl$rhoend, ctrl$rhoend <= ctrl$rhobeg)
    if (verb)
        cat("RHOBEG = ", ctrl$rhobeg,", RHOEND = ", ctrl$rhoend, "\n")

    ## I think this check should now be unnecessary but I'm not sure. DMB
    if (any(rng < 2 * ctrl$rhobeg)) {
        warning("All upper - lower must be >= 2*rhobeg. Changing rhobeg") 
        rhobeg <- 0.2 * min(upper-lower)
    }

    ctrl$wsize <- (ctrl$npt + 5L) * (ctrl$npt + n) + 3L * (n * (n + 5L))/2L
  
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
    if (ctrl$maxfun < 10 * n^2)
        warning("'maxfun' less than 10*length(par)^2 not recommended.")

    ## Pass control as an environment with easier access to name/value pairs
    ctl <- new.env(parent = emptyenv())
    lapply(names(ctrl), function(nm) assign(nm, ctrl[[nm]], envir = ctl))

    out <- .Call("bobyqa_c", par, lower, upper,
               function(par) fn(par, ...), ctl, new.env(),
               PACKAGE = "minqa")
  
    class(out) <- "bobyqa"
    out
}

print.bobyqa <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("bobyqa results\n")

  cat("parameter estimates:", toString(x$par), "\n")
  cat("function evaluations:", toString(x$feval), "\n")
  cat("objective function value:", toString(x$fval), "\n")
  
  invisible(x)
}
