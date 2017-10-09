## @knitr polyobjq

polyobjq <- function(x, penfactor=0, epsilon=0) {
  # negative area + penfactor*(sum(squared violations))
  nv = (length(x)+3)/2 # number of vertices
  area  <-  polyarea(x) # negative area
  f <- -area
  XY <- polypar2XY(x)
  dist2 <- polydistXY(XY)
  viol <- dist2[which(dist2 > 1)] - 1.0
  f <- f + penfactor * sum(viol)
  slacks <- 1.0 + epsilon - dist2 # slack vector
  if (any(slacks <= 0)) { 
    attr(f,"area") <- -area
  } # in case of step into infeasible zone
  else {  
    attr(f,"area") <- area
  }
  attr(f,"minslack") <- min(slacks)
  f
}
