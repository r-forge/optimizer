## @knitr PolyTrack

# 160807 -- need to look at symmetry of result
#           Alternative parametrization

nvex <- 6 # default to hexagon

# Now try to ONLY plot "best so far" polygons


addplot <- function(penv, x, f, area) {
   val <- 1 # OK if returns 0
   nplot <- 5
   ncol <- dim(penv$psave)[2] - 2
#- To add point to the penv$psave matrix and plot the
#- last nplot polygons
   npoint <- dim(penv$psave)[1]
   if (area > penv$besta) {
     penv$besta <- area
     penv$psave <- rbind(penv$psave, c(x, f, area))
     if (npoint > nplot + 1) {
       nset <- nplot
       prh <- penv$psave[1,1:ncol]
       xyrh <- polypar2XY(prh)
       plot.new()
       plot.window(xlim = c(-0.5, 0.7), ylim=c(-0.1, 1.1))
       txt<-paste("Area to reg. polygon = ",area/penv$regarea,sep='')
       title(txt)
       box()
       axis(1)
       axis(2)
       points(xyrh, col='pink', type='l', lwd=2, )
       colrs <- hsv(0.6, (1:nset)/nset, 1, 0.5)
       for (ii in 1:nset) {
         ppoint <- penv$psave[npoint-ii+1, 1:ncol]
         xy <- polypar2XY(ppoint)
         points(xy, col=colrs[ii], type='l', lwd=2)
       }
     }
   }
   val <- 0
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


## @knitr polyobjbig

polyobjbig <- function(x, bignum=1e10, epsilon=0) {
 # Put objective to bignum when constraints violated
 nv = (length(x)+3)/2 # number of vertices
 area <- polyarea(x)
 d2 <- c(x[1:(nv-1)]^2, polypardist2(x)) # distances
 slacks <- 1.0 + epsilon - d2 # slack vector
 if (any(d2 >=1) ) { 
     f <- bignum 
     attr(f,"area") <- -area
 } else { 
    f <-  -area 
    attr(f,"area") <- area
 } # negative area
 attr(f,"minslack") <- min(slacks)
 f
}


## @knitr polyobju

polyobju <- function(x, penfactor=1e-5, epsilon=0, penv) {
  # polyobj with radial parameters constrained by log barrier
  # epsilon <- 0
  bignum <- 1e+20
  # (negative area) + penfactor*(sum(squared violations))
  nv = (length(x)+3)/2 # number of vertices
  area <-  polyarea(x) 
  f <- -area # negative area
  dist2 <- polypardist2(x) # from radial coords, excluding radii (bounded)
  dist2 <- c(x[1:(nv-1)]^2, dist2) # Add in radials. Note the squared distances used
  slacks <- 1.0 + epsilon - dist2 # slack vector
  if (any(slacks <= 0)) { 
    f <- bignum 
    area <- -area # invalid polygon
    attr(f,"area") <- area
  } # in case of step into infeasible zone
  else {  
    f <- f - penfactor*sum(log(slacks)) 
    attr(f,"area") <- area
    addplot(penv, x, f, area)
  }
  attr(f,"minslack") <- min(slacks)
  f
}

## @knitr polygradu

polygradu <- function(x, penfactor=1e-8, epsilon=0, penv) {
  nv <- (length(x)+3)/2
  l8 <- nv - 3 # end of radii params
  # epsilon <- 0
  # (negative area) + penfactor*(sum(squared violations))
  nn <- length(x)
  gg <- rep(0, nn)
  dist2 <- polypardist2(x) # from radial coords, excluding radii (bounded)
  dist2 <- c(x[1:(nv-1)]^2, dist2)
  slacks <- 1.0 + epsilon - dist2 # slack vector
  
  if (any(slacks <= 0)) { # Leave gradient at 0, rely on bignum in polyobju
    cat("polygrad: Infeasible parameters at\n")
    print(x)
    oldw <- getOption("warn")
    options(warn = -1)
    warning("Polygradu -- Infeasible")  
    options(warn = oldw)
  } else { 
    for (ll in 3:nv) {
      ra<-x[ll-1]
      rb<-x[ll-2]
      abangle <- x[l8 + ll]
      # are is 0.5*ra*rb*sin(abangle)
      gg[ll-2] <- gg[ll-2] - 0.5*ra*sin(abangle)
      gg[ll-1] <- gg[ll-1] - 0.5*rb*sin(abangle)
      gg[ll+l8] <- gg[ll+l8] - 0.5*ra*rb*cos(abangle)
    }
  }
  ll <- 0
  # components from radial parameter constraints (upper bounds)
  for (ii in 1:(nv-1)){
    ll <- ll+1
    gg[ll] <- gg[ll] + 2*penfactor*x[ll]/slacks[ll]
  }
  # components from other distances
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


## @knitr polysetup
polysetup <- function(nv, defsize=0.98){
# Function to set up animation of the "largest small polygon"
#   problem. This attempts to find the polygon in nv vertices
#   that has the largest area inside the polygon subject to
#   the constraint that no two vertices are more than 1 unit
#   distant from each other.
# Ref. Graham, "The largest small hexagon" ....???
    cat("polysetup with ",nv," vertices\n")
#    nv <- readline("number of vertices = ")
    nvmax <- 100 # Arbitrary limit -- change ??
    cat("nv, nvmax:",nv, nvmax, "\n")
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
   # defsize <- 1
    regangle <- pi/nv #  pi/no. of vertices
# test to define polygon
    q5<-defsize*sin(regangle) # REM regangle/nv = alpha
    b<-rep(NA,n)
#    x <- rep(NA, nv)
#    y <- rep(NA, nv)
#    x[1] <- 0
#    y[1] <- 0
#    x[2] <- q5
#    y[2] <- 0
    b[1]<-q5
    q1 <- q5
    q2 <- 0 # x2 and y2
    l8 <- nv - 3 # offset for indexing
    for (ll in 3:nv){
        b[ll+l8] <- regangle
        q1 <- q1+q5*cos(2*(ll-2)*regangle)
        q2 <- q2+q5*sin(2*(ll-2)*regangle)
#        x[ll]<-q1
#        y[ll]<-q2
        b[ll-1]<-sqrt(q1*q1+q2*q2)
    }
#    par0 <- b # return the parameters as par0
    res <- list(par0 = b, lb = lb, ub =ub)
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

