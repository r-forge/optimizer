\name{ControlParams}
\alias{ControlParams}
\encoding{UTF-8}
\title{set control defaults}
\concept{minimization}
\concept{maximization}
\description{
  Set control defaults.
}
\usage{
   ControlParams(LBFGS = FALSE, memory = 11, SubCheck = 8, SubSkip = 4,
                 eta0 = 0.001, eta1 = 0.9, eta2 = 1e-10, AWolfe =
                 FALSE, AWolfeFac = 0.001, Qdecay = 0.7, nslow = 1000,
                 StopRule = TRUE, StopFac = 0, PertRule = TRUE, eps =
                 1e-06, egrow = 10, QuadStep = TRUE, QuadCutOff =
                 1e-12, QuadSafe = 1e-10, UseCubic = TRUE, CubicCutOff
                 = 1e-12, SmallCost = 1e-30, debug = FALSE, debugtol =
                 1e-10, step = 0, maxit = .Machine$integer.max, ntries
                 = 50, ExpandSafe = 200, SecantAmp = 1.05, RhoGrow = 2,
                 neps = 5, nshrink = 10, nline = 50, restart_fac = 6,
                 feps = 0, nan_rho = 1.3, nan_decay = 0.1)
}
\arguments{
    \item{LBFGS}{See Hager and Zhang manual.}
    \item{memory}{See Hager and Zhang manual.}
    \item{SubCheck}{See Hager and Zhang manual.}
    \item{SubSkip}{See Hager and Zhang manual.}
    \item{eta0}{See Hager and Zhang manual.}
    \item{eta1}{See Hager and Zhang manual.}
    \item{eta2}{See Hager and Zhang manual.}
    \item{AWolfe}{See Hager and Zhang manual.}
    \item{AWolfeFac}{See Hager and Zhang manual.}
    \item{Qdecay}{See Hager and Zhang manual.}
    \item{nslow}{See Hager and Zhang manual.}
    \item{StopRule}{See Hager and Zhang manual.}
    \item{StopFac}{See Hager and Zhang manual.}
    \item{PertRule}{See Hager and Zhang manual.}
    \item{eps}{See Hager and Zhang manual.}
    \item{egrow}{See Hager and Zhang manual.}
    \item{QuadStep}{See Hager and Zhang manual.}
    \item{QuadCutOff}{See Hager and Zhang manual.}
    \item{QuadSafe}{See Hager and Zhang manual.}
    \item{UseCubic}{See Hager and Zhang manual.}
    \item{CubicCutOff}{See Hager and Zhang manual.}
    \item{SmallCost}{See Hager and Zhang manual.}
    \item{debug}{See Hager and Zhang manual.}
    \item{debugtol}{See Hager and Zhang manual.}
    \item{step}{See Hager and Zhang manual.}
    \item{maxit}{See Hager and Zhang manual.}
    \item{ntries}{See Hager and Zhang manual.}
    \item{ExpandSafe}{See Hager and Zhang manual.}
    \item{SecantAmp}{See Hager and Zhang manual.}
    \item{RhoGrow}{See Hager and Zhang manual.}
    \item{neps}{See Hager and Zhang manual.}
    \item{nshrink}{See Hager and Zhang manual.}
    \item{nline}{See Hager and Zhang manual.}
    \item{restart_fac}{See Hager and Zhang manual.}
    \item{feps}{See Hager and Zhang manual.}
    \item{nan_rho}{See Hager and Zhang manual.}
    \item{nan_decay}{See Hager and Zhang manual.}
}
\value{
  \code{ControlParams} returns the default control settings for optimization tools.
}
\keyword{nonlinear}
\keyword{optimize}

