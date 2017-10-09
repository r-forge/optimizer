## @knitr polysetup

polysetup <- function(nv, defsize=0.98){
# Function to set up animation of the "largest small polygon"
#   problem. This attempts to find the polygon in nv vertices
#   that has the largest area inside the polygon subject to
#   the constraint that no two vertices are more than 1 unit
#   distant from each other.
#   Ref. Graham, "The largest small hexagon" 
#    (J. Combinatorial Theory (A), vol. 18, pp. 165-170, 1975)

    nvmax <- 100 # Arbitrary limit -- change ??

    if (nv > nvmax) { stop("Too many vertices for polygon") }

    mcon <- (nv-2)*(nv-1)/2 # Number of distance constraints

    n <- 2*nv - 3 # Number of parameters in the problem

    # Thus we use a vector b[] of length n
    # Note that we use RADIAL coordinates to simplify the
    # optimization, but convert to cartesian to plot them
    # First point is always at the origin (0,0) cartesian
    # Second point is at (b[1],0) in both cartesian or polar
    # where cartesian is (x, y) and radial is (radius, angle)
    # Choice: angle in radians. ??
    # There are 2*nv cartesian coordinate values
    # i.e., (x, y) for nv point
    # But first point is (0,0) and second has angle 0
    #   since point 2 fixed onto x axis (angular coordinate 0).
    # So b[1] ... b[nv-1] give radial coordinates of points 2:nv
    # and b[nv] ... b[2*nv-3] give angle coordinates of points 3:nv
    # ?? not needed LET L8=nv-3: REM so l+l8 indexes angles as l=3..nv
    # Distances between points can be worked out by cosine rule for
    # triangles i.e. D = sqrt(ra^2 + rb^2 - 2 ra rb cos(angle)
    # Now set lower and upper bounds

    lb <- rep(0, n) # all angles and distances non-negative

    ub <- c(rep(1, (nv-1)), rep(pi, (nv-2))) # distances <=1, angles <= pi

    # if we have angles > pi, then we are reflecting the polygon about an edge
    # set inital parameters to a regular polygon of size .98
#    defsize <- 0.98
    regangle <- pi/nv #  pi/(no. of vertices)
# test to define polygon
    q5 <- defsize * sin(regangle) # REM regangle/nv = alpha
    b <- rep(NA,n)
#    x <- rep(NA, nv)
#    y <- rep(NA, nv)
#    x[1] <- 0
#    y[1] <- 0
#    x[2] <- q5
#    y[2] <- 0
    b[1] <- q5
    q1 <- q5
    q2 <- 0 # x2 and y2
    l8 <- nv - 3 # offset for indexing
    for (ll in 3:nv){
        b[ll+l8] <- regangle
        q1 <- q1 + q5*cos(2 * (ll-2) * regangle)
        q2 <- q2 + q5 * sin(2 * (ll-2) * regangle)
#        x[ll]<-q1
#        y[ll]<-q2
        b[ll-1] <- sqrt(q1*q1 + q2*q2) # radius
    }
    res <- list(par0 = b, lb = lb, ub = ub)
}

## @knitr polypar2XY

polypar2XY <- function(b) {
# converts radial coordinates for polygon into Cartesian coordinates
#  that are more suitable for plotting
    nv <- (length(b)+3)/2
    l8 <- nv - 3 # offset for indexing
    x <- rep(NA, nv+1)
    y <- rep(NA, nv+1)
    # One extra point to draw polygon (return to origin)
    x[1] <- 0
    y[1] <- 0
    x[2] <- b[1]
    y[2] <- 0
    cumangle <- 0 # Cumulative angle of points so far
    q5 <- b[1]
    q1 <- q5 # x2
    q2 <- 0 #  y2
    for (ll in 3:nv){
        cumangle <- cumangle + b[ll+l8]
        cradius <- b[ll-1]
        q1 <- cradius*cos(cumangle)
        q2 <- cradius*sin(cumangle)
        x[ll]<-q1
        y[ll]<-q2
    }
    x[nv+1] <- 0 # to close the polygon
    y[nv+1] <- 0    
    XY <- list(x=x, y=y)
    XY
}

## @knitr polyarea

polyarea<-function(b) {
   # compute area of a polygon defined by radial coordinates
   # This IGNORES constraints
   nv <- (length(b)+3)/2
   area <- 0 
   l8 <- nv-3
   for (l in 3:nv){ # nv - 2 triangles
      q1 <- b[[l-2]] # side 1
      q2 <- b[[l-1]] # side 2
      q3 <- b[[l+l8]] # angle
      atemp <- q1*q2*sin(q3)
      area <- area + atemp
   }
   area <- area * 0.5
   area
}

## @knitr polydistXY

polydistXY <- function(XY) {
#   compute point to point distances from XY data
   nv <- length(XY$x)-1
   ncon <- (nv - 1)*(nv)/2
   dist2 <- rep(NA, ncon) # squared distances   
   ll <- 0 # index of constraint
   for (i in 1:(nv-1)){
      for (j in ((i+1):nv)){
         xi <- XY$x[i]
         xj <- XY$x[j]
         yi <- XY$y[i]
         yj <- XY$y[j]
         dd <- (xi-xj)^2 + (yi-yj)^2
         ll <- ll + 1
         dist2[ll] <- dd
      }
   }        
   dist2
}

## @knitr polypar2distXY

polypar2distXY <- function(pars) {
# compute the pairwise distances using two calls
   nv <- (length(pars) + 3)/2
   XY <- polypar2XY(pars)
   dist2 <- polydistXY(XY)
}


## @knitr polypardist2

polypardist2 <- function(b) {
# compute the pairwise distances for non-radii lines
   nv <- (length(b) + 3)/2 
   l8 <- nv - 3 # end of radii params
   ll <- 0 # count the distances (non-radii ones)
   sqdist <- rep(NA, (nv-1)*(nv-2)/2)
   for (ii in 2:(nv-1)){
      for (jj in (ii+1):nv) {
          ra <- b[ii-1]
          rb <- b[jj-1]
          angleab <- 0
          for (kk in (ii+1):jj) { angleab <- angleab + b[kk+l8] }
          d2 <- ra*ra+rb*rb -2*ra*rb*cos(angleab) # Cosine rule for squared dist
          ll <- ll+1
          sqdist[[ll]] <- d2
      }
   }  
   sqdist
}

## @knitr polyobj

polyobj <- function(x, penfactor=1e-8, epsilon=0) {
# log barrier objective function for small polygon
# epsilon <- 0
 bignum <- 1e+20
 # (negative area) + penfactor*(sum(squared violations))
 nv = (length(x)+3)/2 # number of vertices
 area <- polyarea(x) # area
 f <- - area 
 dist2 <- polypardist2(x) # from radial coords, excluding radii (bounded)
 slacks <- 1.0 + epsilon - dist2 # slack vector
 if (any(slacks <= 0)) { 
#     cat("polygrad: Infeasible parameters at\n")
#     print(x)
     f <- bignum 
     area <- -area # to code for infeasible and avoid plotting
 } # in case of step into infeasible zone
 else {  f <- f - penfactor*sum(log(slacks)) }
 attr(f,"area") <- area
 attr(f,"minslack") <- min(slacks)
 f
}


## @knitr polygrad

polygrad <- function(x, penfactor=1e-8, epsilon=0) {
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
    stop("polygrad: Infeasible") 
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


## @knitr polyobjp

polyobjp <- function(x, penfactor=1e-8, epsilon=0) {
# log barrier objective function for small polygon
# epsilon <- 0
 bignum <- 1e+20
 # (negative area) + penfactor*(sum(squared violations))
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

## @knitr PolyTrack

library(R6)
library(TeachingDemos)

 nvex <- 6 # default value -- need to check if we can change 

PolyTrack <- R6Class("PolyTrack",
  public = list(
    parms = list(),
    maxviol = list(),
    areas = list(),
    fvals =list(),
    nv = nvex,
    PlotIt = TRUE,
    Delay = 0.25,
    nPolys = 5,
    add = function(p,v,a, fval) { # add points of polygon and area
      i <- length(self$parms) + 1
      self$parms[[i]] <- p # the points
      self$maxviol[[i]] <- v # maximum violation
      self$areas[[i]] <- a # the area
      self$fvals[[i]] <- fval # objective
      if(self$PlotIt) { # here PlotIt in environment is TRUE so we'll likely always do this
        self$PlotPolys() # plot all polygons to date, then wait
        Sys.sleep(self$Delay)
      }
      return(a)
    },
    PlotPolys = function(i=-1) { # to draw the polygons
      if(i<0) i <- length(self$parms)
      if(i==0) return()
      cols <- hsv(0.6, (1:self$nPolys)/self$nPolys, 1)
      # sets up a vector of colours. In this case we want gradual fade-out
      # of the older polygons so we can see the latest the best
      start <- pmax(1, i-self$nPolys+1)
      plotParms <- self$parms[seq(start,i)]
      n <- length(plotParms)
      if(n < self$nPolys) cols <- tail(cols, n)
      coords <- lapply(plotParms, function(x) polypar2XY(x))
      plot.new()
      plot.window( xlim=do.call(range, lapply(coords, function(xy) xy$x)),
                   ylim=do.call(range, lapply(coords, function(xy) xy$y)),
                   asp=1)
      for(ii in seq_len(n)) { # draw the edges of polygons in the set
        polygon(coords[[ii]]$x, coords[[ii]]$y, border=cols[ii], lwd=3)
      }
      # Here add display of area found
      carea <- max(unlist(self$areas))
#      txt <- paste("Max polygon area =",carea,"  last=",unlist(self$areas)[-1])
      txt <- paste("Max polygon area =",carea)
      title(main=txt)
      title(sub=paste("Max violation=",self$maxviol[[i]],"  obj.fn.=",self$fvals[[i]]))
    }
  ) # end of public list, no private list
# NOTE: Need to change nv to whatever is current value
                     
)

