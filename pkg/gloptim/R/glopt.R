##
##  g l o p t . R
##


glopt <- function(fn, lb, ub,
                  x0 = NULL, rand = FALSE,
                  maxiter = NULL, popsize = NULL,
                  g = NULL, gr = NULL,
                  incl = NULL, excl = NULL) {
    all_methods <-
        c("deoptim", "cppdeoptim", "deoptimr",          # **DE**
          "deopt", "simplede", "simpleea",              # **EA**
          "gensa", "ga", "genoud",                      # **GA**
          "pso", "psopt", "hydropso",                   # **PSO**
          "direct", "crs2lm", "isres", "stogo",         # **NLoptr**
          "cmaoptim", "cmaes", "cmaesr", "purecmaes",   # **CMA-ES**
          "malschains",  "ceimopt",                     # **CE**
          "smco", "soma")                               # --others--

    if (is.null(incl))  incl <- all_methods
    if (!is.null(excl)) incl <- setdiff(incl, excl)

    mthd <- c(); fminimum <- c(); thetime = c()
    for (m in incl) {
        tm <- system.time(
            sol <- gloptim(fn, lb, ub, x0 = x0, method=m, gr=gr)
        )
        mthd <- c(mthd, m); fminimum <- c(fminimum, sol$fmin)
        thetime <- c(thetime, unname(tm["elapsed"]))
    }
    # Sort according to minimal value and return data frame
    o <- order(fminimum, na.last = TRUE)
    data.frame(method=mthd[o], fmin=fminimum[o], time=thetime[o])
}
