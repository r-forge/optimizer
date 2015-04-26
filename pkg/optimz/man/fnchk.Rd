\name{fnchk}
\alias{fnchk}
\title{Run tests, where possible, on user objective function}
\description{
   \code{fnchk} checks a user-provided R function, \code{ffn}. 
}
\usage{
   fnchk(xpar, ffn, trace=0, ... )
}
\arguments{
     
    \item{xpar}{
        the (double) vector of parameters to the objective funcion
        }
    \item{ffn}{
        a user-provided function to compute the objective function
        }
    \item{trace}{
	set >0 to provide output from fnchk to the console, 0 otherwise
        }
    \item{\dots}{
        optional arguments passed to the objective function.
        }
}
\details{
\code{fnchk} attempts to discover various errors in function setup in user-supplied
functions primarily intended for use in optimization calculations. There are always 
more conditions that could be tested!    
}
\value{
  The output is a list consisting of 
list(fval=fval, infeasible=infeasible, excode=excode, msg=msg)
  \item{fval}{The calculated value of the function at parameters \code{xpar} if the function
      can be evaluated.}
  \item{infeasible}{FALSE if the function can be evaluated, TRUE if not.}
  \item{excode}{An exit code, which has a relationship to }
  \item{msg}{A text string giving information about the result of the function check: Messages and
    the corresponding values of \code{excode} are:
    \itemize{
        \item{fnchk OK;}{ \code{excode} = 0; 
                       \code{infeasible} = FALSE}
        \item{Function returns INADMISSIBLE;}
                    { \code{excode} = -1; \code{infeasible} = TRUE}
        \item{Function returns a vector not a scalar;}
                    { \code{excode} = -4; \code{infeasible} = TRUE}
        \item{Function returns a list not a scalar;}
                    { \code{excode} = -4; \code{infeasible} = TRUE}
        \item{Function returns a matrix list not a scalar;}
                    { \code{excode} = -4; \code{infeasible} = TRUE}
        \item{Function returns an array not a scalar;}
                    { \code{excode} = -4; \code{infeasible} = TRUE}
        \item{Function returned not length 1, despite not vector, matrix or array;}
                    { \code{excode} = -4; \code{infeasible} = TRUE}
        \item{Function returned non-numeric value; \code{excode} = 0;}
                    { \code{excode} = -1; \code{infeasible} = TRUE}
        \item{Function returned Inf or NA (non-computable);}
                    { \code{excode} = -1; \code{infeasible} = TRUE}
    }
  }
}


\author{
    John C. Nash <nashjc@uottawa.ca>
}

\examples{
cat("all examples in tests/ subdirectory\n")

}
\keyword{optimize}

