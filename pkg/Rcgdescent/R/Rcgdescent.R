## The wrapped C code sets a large number of control parameters in a complicated
## structure.  The ControlParams function allows the user to set these
## parameters from R, if desired.  The defaults are often fine.
ControlParams <- function(LBFGS = FALSE,
  memory = 11,
  SubCheck = 8,
  SubSkip = 4,
  eta0 = 0.001,
  eta1 = 0.900,
  eta2 = 1.0e-10,
  AWolfe = FALSE,
  AWolfeFac = 1.0e-3,
  Qdecay = .7,
  nslow = 1000,
  StopRule = TRUE,
  StopFac = 0.e-12,
  PertRule = TRUE,
  eps = 1.0e-6,
  egrow = 10,
  QuadStep = TRUE,
  QuadCutOff = 1.0e-12,
  QuadSafe = 1e-10,
  UseCubic = TRUE,
  CubicCutOff = 1e-12,
  SmallCost = 1e-30,
  debug = FALSE,
  debugtol = 1e-10,
  step = 0.000,
  maxit = .Machine$integer.max,
  ntries = 50,
  ExpandSafe = 200,
  SecantAmp = 1.05,
  RhoGrow = 2.0,
  neps = 5,
  nshrink = 10,
  nline = 50,
  restart_fac = 6,
  feps = 0.000,
  nan_rho = 1.3,
  nan_decay = 0.1) {

  params <- list(LBFGS = as.logical(LBFGS),
    memory = as.integer(memory),
    SubCheck = as.integer(SubCheck),
    SubSkip = as.integer(SubSkip),
    eta0 = as.double(eta0),
    eta1 = as.double(eta1),
    eta2 = as.double(eta2),
    AWolfe = as.integer(AWolfe),
    AWolfeFac = as.double(AWolfeFac),
    Qdecay = as.double(Qdecay),
    nslow = as.integer(nslow),
    StopRule = as.logical(StopRule),
    StopFac = as.double(StopFac),
    PertRule = as.logical(PertRule),
    eps = as.double(eps),
    egrow = as.double(egrow),
    QuadStep = as.logical(QuadStep),
    QuadCutOff = as.double(QuadCutOff),
    QuadSafe = as.double(QuadSafe),
    UseCubic = as.logical(UseCubic),
    CubicCutOff = as.double(CubicCutOff),
    SmallCost = as.double(SmallCost),
    debug = as.logical(debug),
    debugtol = as.double(debugtol),
    step = as.double(step),
    maxit = as.integer(maxit),
    ntries = as.integer(ntries),
    ExpandSafe = as.double(ExpandSafe),
    SecantAmp = as.double(SecantAmp),
    RhoGrow = as.double(RhoGrow),
    neps = as.integer(neps),
    nshrink = as.integer(nshrink),
    nline = as.integer(nline),
    restart_fac = as.double(restart_fac),
    nan_rho = as.double(nan_rho),
    nan_decay = as.double(nan_decay))

  for (prm in params) {
    stopifnot(length(prm) == 1)
  }
  class(params) <- "ControlParams"
  return(params)
}

FunctionWrapper <- function(f, ...) {
  wrapped <- function(x) {
    return(f(x, ...))
  }
  ans <- list(function.name = 'wrapped', env  = environment(wrapped)) 
}

Rcgdescent <- function(par,
                        fn,
                        gr,
                        control = list(maxit=500),
                        ...) {
  ## Args:
  ##   par: A numeric vector containing the values of 'x' where the
  ##     optimization should begin.
  ##   fn: A function to be optimized.  The return value should be a
  ##     numeric scalar.  The first argument is the target of the optimization
  ##     routine.
  ##   gr: A function returning the gradient of 'fn'.  The
  ##     arguments of 'gr' must match 'fn' and be in the same
  ##     order.
  ##   control: A list of control parameters for fine tuning the optimization.
  ##     In most cases the default values are fine.  See help(ControlParams) for
  ##     details.
  ##   ...: If fn and gr take additional parameters, they can be
  ##     supplied by name here.
  stopifnot(is.function(fn))
  stopifnot(is.function(gr))

  ctrl <- ControlParams()
  stopifnot(inherits(ctrl, "ControlParams"))
  ncontrol <- names(control)
  nctrl <- names(ctrl)
  for (onename in ncontrol) {
     if (onename %in% nctrl) {
       if (! is.null(control[onename]) || ! is.na(control[onename]) )
       ctrl[onename]<-control[onename]
     }
  }

  target <- FunctionWrapper(fn, ...)
  grad <- FunctionWrapper(gr, ...)
  cans <- .Call("cgminu_wrapper", target, grad, par, ctrl)
  names(cans) <- c("par", "f", "eval.f", "eval.grad", "message", "converged")
  cgdmsg <- c("convergence tolerance satisfied", 
              "change in func <= feps*|f|", 
              "total number of iterations exceeded maxit", 
              "slope always negative in line search", 
              "number of line search iterations exceeds nline", 
              "search direction not a descent direction", 
              "excessive updating of eps", 
              "Wolfe conditions never satisfied", 
              "debugger is on and the function value increases", 
              "no cost or gradient improvement in 2n + Parm->nslow iterations", 
              "out of memory", 
              "function nan or +-INF and could not be repaired", 
              "invalid choice for memory parameter")
  cans$message <- cgdmsg[(cans$converged + 1)]

  ans <- list(par=cans$par, value=cans$f, counts=c(cans$eval.f, cans$eval.grad),
                convergence=cans$converged, message=cans$message)
  return(ans)
}
                         
