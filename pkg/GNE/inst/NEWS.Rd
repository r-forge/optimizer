\name{NEWS}
\title{\pkg{GNE} News}
\encoding{UTF-8}


\section{CHANGES IN \pkg{GNE} VERSION 0.99-5}{
  \subsection{NEW}{
    \itemize{
      \item{remove the dependency on sensitivity package by ignoring corresponding test files.}
      \item{use Rd-type NEWS file.}
      \item{change vignette format from Rnw to Rmd.}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item{remove a C warning ("this 'if' clause does not guard... [-Wmisleading-indentation]") in SSR.c.}
      \item{remove an error ("DLL requires the use of native symbols") for C entry points.}
      \item{remove a warning ("Recycling array of length 1 in vector-array arithmetic is deprecated.").}
    }
  }
  
}

\section{CHANGES IN \pkg{GNE} VERSION 0.99-4}{
  \subsection{BUG FIXES}{
    \itemize{
      \item{remove a Warning on the vignette 'GNE-howto.pdf'.}
      \item{remove an ERROR on bad use of "class()==".}
    }
  }
}

\section{CHANGES IN \pkg{GNE} VERSION 0.99-3}{
  \subsection{BUG FIX}{
    \itemize{
      \item{remove an ERROR on bad use of "class()==".}
    }
  }
}

\section{CHANGES IN \pkg{GNE} VERSION 0.99-2}{
  \subsection{NOTE FIXES}{
    \itemize{
      \item{remove a NOTE on dynamic symbol.}
      \item{remove a NOTE when rebuilding vignette.}
    }
  }
}

\section{CHANGES IN \pkg{GNE} VERSION 0.99-2}{
  \subsection{NOTE FIXES}{
    \itemize{
      \item{remove a NOTE on dynamic symbol.}
      \item{remove a NOTE when rebuilding vignette.}
    }
  }
}

\section{CHANGES IN \pkg{GNE} VERSION 0.99-1}{
  \subsection{BUG FIXES}{
    \itemize{
      \item{update for R 3.3.0.}
      \item{change title and remove a latex error.}
    }
  }
}

\section{CHANGES IN \pkg{GNE} VERSION 0.99-0}{
  \subsection{BUG FIXES}{
    \itemize{
      \item{remove eqsolve manual.}
      \item{remove C-level warnings on Mac OS Mavericks.}
      \item{remove ampersands.}
    }
  }
}

\section{START \pkg{GNE} VERSION 0.98}{
  \subsection{NEW}{
    \itemize{
      \item{first release on CRAN after removing unnecessary file/directory.}
    }
  }
}

\section{DEVEL \pkg{GNE} VERSION 0.97}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item{add new algorithms and change the arguments of functions GNE.nseq and GNE.ceq.}
      \item{update vignette, manuals and tests files.}
      \item{depends now the SQUAREM package.}
      \item{update the benchmark functions for GNE.nseq.}
      \item{add new constrained equation algorithms.}
      \item{start Stackelberg equilibrium (experimental stage).}
    }
  }
}

\section{DEVEL \pkg{GNE} VERSION 0.96}{
  \subsection{DEFUNC}{
    \itemize{
      \item{remove PROBj2k***, GrPROBj2k***, GrGrPROBj2k*** functions.}
    }
  }
}

\section{DEVEL \pkg{GNE} VERSION 0.95}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item{Wide and deep rebuild of the package. Now use
the BB, alabama, nleqslv packages. }
    }
  }
  \subsection{DEFUNC}{
    \itemize{
      \item{deleted functions are PROBj2k, GrPROBj2k, GrGrPROBj2k and
depreciated functions are NewtonKKT, fixedpoint, minGap.}
    }
  }
}

\section{DEVEL \pkg{GNE} VERSION 0.9}{
  \subsection{NEW}{
    \itemize{
      \item{Initial release on R-forge.}
    }
  }
}

