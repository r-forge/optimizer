## @knitr polyobjbig

polyobjbig <- function(x, bignum=1e10, epsilon=0) {
  # Put objective to bignum when constraints violated
  nv = (length(x)+3)/2 # number of vertices
  area <- polyarea(x)
  d2 <- c(x[1:(nv-1)]^2, polypardist2(x)) # distances
  slacks <- 1.0 + epsilon - d2 # slack vector
  if (any(d2 >=1) ) { 
    f <- bignum 
    attr(f,"area") <- -area # set negative so it is not recorded
  } else { 
    f <-  -area 
    attr(f,"area") <- area
  } # negative area
  attr(f,"minslack") <- min(slacks)
  f
}
