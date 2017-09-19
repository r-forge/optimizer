## @knitr polyobjp

polyobjp <- function(x, penfactor=1e-8, epsilon=0) {
  # log barrier objective function for small polygon
  # epsilon <- 0
  bignum <- 1e+20
  # (negative area) + penfactor*(sum(log(slacks)))
  nv = (length(x)+3)/2 # number of vertices
  cat("x:")
  print(x)
  area <- polyarea(x) # area
  f <- - area 
  dist2 <- polypardist2(x) # from radial coords, excluding radii (bounded)
  slacks <- 1.0 + epsilon - dist2 # slack vector
  if (any(slacks <= 0)) { 
    cat("polyobjp: Infeasible parameters at\n")
    print(x)
    f <- bignum 
    area <- -area # to code for infeasible and avoid plotting
  } # in case of step into infeasible zone
  else {  f <- f - penfactor*sum(log(slacks)) }
  attr(f,"area") <- area
  attr(f,"minslack") <- min(slacks)
  f
}

## @knitr polygradp

polygradp <- function(x, penfactor=1e-8, epsilon=0) {
  # log barrier gradient function for small polygon
  nv <- (length(x)+3)/2
  l8 <- nv - 3 # end of radii params
  # epsilon <- 0
  bignum <- 1e+20
  # (negative area) + penfactor*(sum(squared violations))
  nn <- length(x)
  gg <- rep(0, nn)
  dist2 <- polypardist2(x) # from radial coords, excluding radii (bounded)
  slacks <- 1.0 + epsilon - dist2 # slack vector
  if (any(slacks <= 0)) { 
    cat("polygrad: Infeasible parameters at\n")
    print(x)
    warning("polygrad: Infeasible") 
    return(gg)
  } 
  for (ll in 3:nv) {
    ra<-x[ll-1]
    rb<-x[ll-2]
    abangle <- x[l8 + ll]
    # are is 0.5*ra*rb*sin(abangle)
    gg[ll-2] <- gg[ll-2] - 0.5*ra*sin(abangle)
    gg[ll-1] <- gg[ll-1] - 0.5*rb*sin(abangle)
    gg[ll+l8] <- gg[ll+l8] - 0.5*ra*rb*cos(abangle)
  }
  ll <- 0
  for (ii in 2:(nv-1)){
    for (jj in (ii+1):nv) {
      ll <- ll+1
      ra <- x[ii-1]
      rb <- x[jj-1]
      angleab <- 0
      for (kk in (ii+1):jj) { angleab <- angleab + x[kk+l8] }
      gg[ii-1] <- gg[ii-1] + 2*penfactor*(ra-rb*cos(angleab))/slacks[ll]
      gg[jj-1] <- gg[jj-1] + 2*penfactor*(rb-ra*cos(angleab))/slacks[ll]
      for (kk in (ii+1):jj){
        gg[kk+l8]<-gg[kk+l8]+2*penfactor*ra*rb*sin(angleab)/slacks[ll]
      }
    }
  }
  gg
}
