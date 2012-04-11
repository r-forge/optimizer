\name{optstart}
\alias{optstart}
\title{Create scratchpad for holding optimization control variables.
}
\description{Creates a scratchpad environment OPCON for holding 
optimization control variables. This has .Globalenv as parent
environment.
}
\usage{
optstart(npar)
}
\arguments{
 \item{npar}{The number of parameters to be optimized.}
}
\details{
None
}
\value{
   The name of the environment that is the scratchpad is returned.
}
\author{
John C Nash <nashjc@uottawa.ca>

Maintainer: John C Nash <nashjc@uottawa.ca>
}

\keyword{ optimization }
\seealso{
optimx
}

