\name{nltestprob}
\alias{Chwirut1}
\alias{Chwirut1.f}
\alias{Chwirut1.jac}
\alias{Chwirut1.g}
\alias{Chwirut1.setup}
\alias{Chwirut1.test}
\alias{Chwirut1.res}
\alias{Chwirut1.h}
\alias{DanielWood}
\alias{DanielWood.f}
\alias{DanielWood.jac}
\alias{DanielWood.g}
\alias{DanielWood.setup}
\alias{DanielWood.test}
\alias{DanielWood.res}
\alias{DanielWood.h}
\alias{Ratkowsky2}
\alias{Ratkowsky2.f}
\alias{Ratkowsky2.jac}
\alias{Ratkowsky2.g}
\alias{Ratkowsky2.setup}
\alias{Ratkowsky2.test}
\alias{Ratkowsky2.res}
\alias{Ratkowsky2.h}
\title{NIST nonlinear regression problems recast as function minimization}
\description{
   The NIST nonlinear regression test problems from 

   http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml

   are provided as data sets plus functions to perform the following calculations:

   \describe{
   \item{setup}{ in a function with form NAME.setup}
   \item{test}{ in a function with form NAME.test}
   \item{residuals}{ in a function with form NAME.res}
   \item{sum of squares function}{ in a function with form NAME.f}
   \item{gradient of the sum of squares}{ in a function with form NAME.g}
   \item{Jacobian of the residuals}{ in a function with form NAME.jac}
   \item{data}{ in a data frame with a given NAME}

   }

   The problems are as follows:

  \describe{
  \item{\code{Chwirut1}: Ultrasonic calibration study 1}{ 
    The \code{Chwirut1} data frame has 214 rows and 2 columns.

  \describe{
    \item{y}{
      A numeric vector of ultrasonic response values
    }
    \item{x}{
      A numeric vector or metal distance values
    }
  }

  These data are the result of a NIST study involving
  ultrasonic calibration.  The response variable is
  ultrasonic response, and the predictor variable is
  metal distance.

  Source:
  Chwirut, D., NIST (197?).   Ultrasonic Reference Block Study. 
  } 

  \item{ \code{DanielWood}: Radiated energy}{
The \code{DanielWood} data frame has 6 rows and 2 columns giving the
energy radiated from a carbon filament versus the absolute temperature
of the filament.

  \describe{
    \item{y}{
      A numeric vector of the energy radiated from a carbon filament
      lamp.
    }
    \item{x}{
      A numeric vector of the temperature of the filament (1000 K).
    }
  }

  These data and model are described in Daniel and Wood
    (1980), and originally published in E.S.Keeping, 
    "Introduction to Statistical Inference," Van Nostrand
    Company, Princeton, NJ, 1962, p. 354.  The response
    variable is energy radiated from a carbon filament
    lamp per cm**2 per second, and the predictor variable
    is the absolute temperature of the filament in 1000
    degrees Kelvin.

Source:
Daniel, C. and F. S. Wood (1980).
Fitting Equations to Data, Second Edition. 
New York, NY:  John Wiley and Sons, pp. 428-431.
}

\item{\code{Ratkowsky2}: Pasture yield data}{
  The \code{Ratkowsky2} data frame has 9 rows and 2 columns.


  \describe{
    \item{y}{
      A numeric vector of pasture yields.
    }
    \item{x}{
      A numeric vector of growing times.
    }
  }


  This model and data are an example of fitting
  sigmoidal growth curves taken from Ratkowsky (1983).
  The response variable is pasture yield, and the
  predictor variable is growing time.

  Source:
Ratkowsky, D.A. (1983).  
Nonlinear Regression Modeling.
New York, NY:  Marcel Dekker, pp. 61 and 88.
}
 
}
}
%%% End of describe
\keyword{datasets}
\keyword{function minimization}
\keyword{optimization}
\keyword{nonlinear regression}
\keyword{nonlinear least squares}

