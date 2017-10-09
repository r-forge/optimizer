## @knitr XRosenbrock.prb
# This is file XRosenbrock.prb
probname <- "XRosenbrock"
probdesc <- " This is the second multidimensional variation of the Rosenbrock function
from Wikipedia's article https://en.wikipedia.org/wiki/Rosenbrock_function. The original
2-parameter function is from Rosenbrock, H.H. (1960). 'An automatic method for finding the
greatest or least value of a function'. The Computer Journal. 3: 175–184. 
doi:10.1093/comjnl/3.3.175. ISSN 0010-4620.
"

#- Note: environment / list "counter" must already exist

if (! exists("counters")) {stop("The 'counters' environment must exist before problems are run.")}

#- setup

XRosenbrock.start <- function(indx) {
  #- indx is character string to allow for more general forms e.g., XRosenbrock
  #- split the index  
  ivec <- strsplit(indx, ":")[[1]] # only want first element of list
  npar <- as.numeric(ivec[1])
  if ( (npar ==2) && (ivec[2]=="HHR")) { start = c(-1.2, 1) } # Traditional start n=2
  else { if (ivec[2] == "pi") { 
           start <- rep(pi, npar) 
         } else { start <- rep(as.numeric(ivec[2]), npar) }
  }
  pstring <- c()
  for (k in 1:npar) {pstring[k] <- paste("p",k,sep='')} # name the parameters
  names(start)<-pstring 
  start
}
XRosenbrock.f <- function (x) 
{
    n <- length(x)
    x1 <- x[2:n]
    x2 <- x[1:(n - 1)]
    sum(100 * (x1 - x2^2)^2 + (1 - x2)^2)
}
XRosenbrock.g <- function (x) 
{
    n <- length(x)
    g <- rep(NA, n)
    g[1] <- 2 * (x[1] - 1) + 400 * x[1] * (x[1]^2 - x[2])
    if (n > 2) {
        ii <- 2:(n - 1)
        g[ii] <- 2 * (x[ii] - 1) + 400 * x[ii] * (x[ii]^2 - x[ii + 
            1]) + 200 * (x[ii] - x[ii - 1]^2)
    }
    g[n] <- 200 * (x[n] - x[n - 1]^2)
    g
}
