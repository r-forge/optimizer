%%% $Id: ENSO.Rd,v 1.3 2003/07/22 19:42:20 bates Exp $
\name{ENSO_opt}
\alias{ENSO}
\alias{ENSO.f}
\alias{ENSO.jac}
\alias{ENSO.g}
\alias{ENSO.setup}
\alias{ENSO.test}
\alias{ENSO.res}
\alias{ENSO.h}
%%% \non_function{}
\title{Atmospheric pressure differences}
\description{
The \code{ENSO} data frame has 168 rows and 2 columns giving atmospheric
pressure differences over time.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{
      A numeric vector of monthly averaged atmospheric pressure
      differences between Easter Island and Darwin, Australia.
    }
    \item{x}{
      A numeric vector of time values.
    }
  }
}
\details{
    The data are monthly averaged atmospheric pressure 
    differences between Easter Island and Darwin, 
    Australia.  This difference drives the trade winds in 
    the southern hemisphere.  Fourier analysis of the data
    reveals 3 significant cycles.  The annual cycle is the
    strongest, but cycles with periods of approximately 44
    and 26 months are also present.  These cycles
    correspond to the El Nino and the Southern Oscillation.
    Arguments to the SIN and COS functions are in radians.
}
\source{
Kahaner, D., C. Moler, and S. Nash, (1989). 
Numerical Methods and Software.  
Englewood Cliffs, NJ: Prentice Hall, pp. 441-445.
}
\examples{
Try <- function(expr) if (!inherits(val <- try(expr), "try-error")) val
plot(y ~ x, data = ENSO)
plot(y ~ x, data = ENSO, type = "l")  # to see the pattern more clearly
Try(fm1 <- nls(y ~ b1 + b2*cos( 2*pi*x/12 ) + b3*sin( 2*pi*x/12 ) 
                   + b5*cos( 2*pi*x/b4 ) + b6*sin( 2*pi*x/b4 )
                   + b8*cos( 2*pi*x/b7 ) + b9*sin( 2*pi*x/b7 ),
            data = ENSO, trace = TRUE,
            start = c(b1 = 11.0, b2 = 3.0, b3 = 0.5, b4 = 40.0, b5 = -0.7,
                      b6 = -1.3, b7 = 25.0, b8 = -0.3, b9 = 1.4)))
Try(fm1a <- nls(y ~ b1 + b2*cos( 2*pi*x/12 ) + b3*sin( 2*pi*x/12 ) 
                    + b5*cos( 2*pi*x/b4 ) + b6*sin( 2*pi*x/b4 )
                    + b8*cos( 2*pi*x/b7 ) + b9*sin( 2*pi*x/b7 ),
             data = ENSO, trace = TRUE, alg = "port",
             start = c(b1 = 11.0, b2 = 3.0, b3 = 0.5, b4 = 40.0, b5 = -0.7,
                       b6 = -1.3, b7 = 25.0, b8 = -0.3, b9 = 1.4)))
Try(fm2 <- nls(y ~ b1 + b2*cos( 2*pi*x/12 ) + b3*sin( 2*pi*x/12 ) 
                   + b5*cos( 2*pi*x/b4 ) + b6*sin( 2*pi*x/b4 )
                   + b8*cos( 2*pi*x/b7 ) + b9*sin( 2*pi*x/b7 ),
            data = ENSO, trace = TRUE,
            start = c(b1 = 10.0, b2 =  3.0, b3 =  0.5, b4 = 44.0, b5 = -1.5,
                     b6 =  0.5, b7 = 26.0, b8 = -0.1, b9 =  1.5)))
Try(fm2a <- nls(y ~ b1 + b2*cos( 2*pi*x/12 ) + b3*sin( 2*pi*x/12 ) 
                    + b5*cos( 2*pi*x/b4 ) + b6*sin( 2*pi*x/b4 )
                    + b8*cos( 2*pi*x/b7 ) + b9*sin( 2*pi*x/b7 ),
             data = ENSO, trace = TRUE, alg = "port",
             start = c(b1 = 10.0, b2 =  3.0, b3 =  0.5, b4 = 44.0, b5 = -1.5,
                     b6 =  0.5, b7 = 26.0, b8 = -0.1, b9 =  1.5)))
Try(fm3 <- nls(y ~ cbind(1, cos( 2*pi*x/12 ), sin( 2*pi*x/12 ), cos( 2*pi*x/b4 ),
                      sin( 2*pi*x/b4 ), cos( 2*pi*x/b7 ), sin( 2*pi*x/b7 )),
            data = ENSO, trace = TRUE,
            start = c(b4 = 40.0, b7 = 25.0), algorithm = "plinear"))
Try(fm4 <- nls(y ~ cbind(1, cos( 2*pi*x/12 ), sin( 2*pi*x/12 ), cos( 2*pi*x/b4 ),
                    sin( 2*pi*x/b4 ), cos( 2*pi*x/b7 ), sin( 2*pi*x/b7 )),
            data = ENSO, trace = TRUE,
            start = c(b4 = 44.0, b7 = 26.0), algorithm = "plinear"))
}
\keyword{datasets}
