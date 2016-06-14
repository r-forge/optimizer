## @knitr polysetup

polysetup <- function(nv, defsize=0.98){
# Function to set up animation of the "largest small polygon"
#   problem. This attempts to find the polygon in nv vertices
#   that has the largest area inside the polygon subject to
#   the constraint that no two vertices are more than 1 unit
#   distant from each other.
# Ref. Graham, "The largest small hexagon" ....???
#    nv <- readline("number of vertices = ")
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
    nv <- (length(b)+3)/2
    l8 <- nv - 3 # offset for indexing
    x <- rep(NA, nv+1)
    y <- rep(NA, nv+1)
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
   nv <- (length(b)+3)/2
   area <- 0 
   l8 <- nv-3
   for (l in 3:nv){
      q1 <- b[[l-2]]
      q2 <- b[[l-1]]
      q3 <- b[[l+l8]]
      atemp <- q1*q2*sin(q3)
      area <- area + atemp
   }
   area <- area * 0.5
   area
}

## @knitr polydistXY

polydistXY <- function(XY) {
#   compute point to point distances from XY data
   nv <- dim(XY)[1]
   ncon <- (nv - 1)*(nv - 2)/2
   dist2 <- rep(NA, ncon) # squared distances   
   ll <- 0 # index of constraint
   for (i in 2:nv){
      for (j in (1:(i-1))){
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
   nv <- (length(pars) + 3)/2
   XY <- polypar2XY(nv, pars)
   dist2 <- polydistXY(nv, XY)
}


## @knitr polypardist2

polypardist2 <- function(b) {
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
# epsilon <- 0
 bignum <- 1e+20
 # 2 * (negative area) + penfactor*(sum(squared violations))
 nv = (length(x)+3)/2 # number of vertices
 f <-  -2 * polyarea(x) # negative area
 dist2 <- polypardist2(x) # from radial coords, excluding radii (bounded)
 slacks <- 1.0 + epsilon - dist2 # slack vector
 if (any(slacks <= 0)) { f <- bignum } 
 else {  f <- f - penfactor*sum(log(slacks)) }
 f
}


## @knitr polygrad

polygrad <- function(x, penfactor=1e-8, epsilon=0) {
 nv <- (length(x)+3)/2
 l8 <- nv - 3 # end of radii params
# epsilon <- 0
 bignum <- 1e+20
 # 2 * (negative area) + penfactor*(sum(squared violations))
 nn <- length(x)
 gg <- rep(0, nn)
 dist2 <- polypardist2(x) # from radial coords, excluding radii (bounded)
 slacks <- 1.0 + epsilon - dist2 # slack vector
 if (any(slacks <= 0)) { stop("Infeasible") } 
 for (ll in 3:nv) {
    ra<-x[ll-1]
    rb<-x[ll-2]
    abangle <- x[l8 + ll]
    # are is ra*rb*sin(abangle)
    gg[ll-2] <- gg[ll-2] - ra*sin(abangle)
    gg[ll-1] <- gg[ll-1] - rb*sin(abangle)
    gg[ll+l8] <- gg[ll+l8] - ra*rb*cos(abangle)
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



## @knitr polyobj2

polyobj2 <- function(x, penfactor=0) {
 epsilon <- 1e-6
 bignum <- 1e+20
 # negative area - penfactor*(sum(log(slacks)))
# but bail and set objective large if any slack <=0
 nv = (length(x)+3)/2 # number of vertices
 f <-  -polyarea(nv, x) # negative area
 XY <- polypar2XY(nv, x)
 dist2 <- polydistXY(nv, XY)
 if (any(dist2 <= 0)) {f <- bignum}
 else { 
    slacks <- 1.0 + epsilon - dist2
    f <- f - penfactor * sum(log(min(slacks)))
 }
 f
}

## @knitr polyobj2a

polyobj2a <- function(x, penfactor=0) {
 # version to allow general unconstrained minimizer
 # However, gradient approximation may cause trouble!
 epsilon <- 1e-9
 bignum <- 1e+20
 # negative area - penfactor*(sum(log(slacks)))
# but bail and set objective large if any slack <=0
 nv <- (length(x)+3)/2 # number of vertices
 f <-  -polyarea(nv, x) # negative area
 XY <- polypar2XY(nv, x)
 dist2 <- polydistXY(nv, XY)
 dist2 <- c(dist2, (x[1:(nv-1)])^2)
 if (any(dist2 >= 1)) {
    if(any(dist2 > 1)) {
      f <- bignum
      cat("VIOLATION:")
      print(dist2)
    } else { # do nothing 
      cat("Max dist\n")
    }
 } else { 
    slacks <- 1.0 + epsilon - dist2
    f <- f - penfactor * sum(log(min(slacks)))
 }
 f
}

## @knitr polyobj3

polyobj3 <- function(x, penfactor=0) {
 # version for DEoptim and similar methods
 epsilon <- 0
 bignum <- 1e+20
 # negative area - penfactor*(sum(log(slacks)))
# but bail and set objective large if any slack <=0
 nv <- (length(x)+3)/2 # number of vertices
 f <-  -polyarea(nv, x) # negative area
 XY <- polypar2XY(nv, x)
 dist2 <- polydistXY(nv, XY)
 dist2 <- c(dist2, (x[1:(nv-1)])^2)
 if (any(dist2 > 1)) {
      f <- bignum
 }
 f
}

## @knitr polyex0

# Example code
nv <- 6
cat("Polygon data:\n")
myhex <- polysetup(nv)
print(myhex)
cat("Area:\n")
myhexa <- polyarea(myhex$par0)
print(myhexa)
cat("XY coordinates\n")
myheXY <- polypar2XY(myhex$par0)
print(myheXY)
plot(myheXY$x, myheXY$y, type="l")
cat("Constraints:\n")
myhexc<-polydistXY(myheXY)
print(myhexc)
cat("Vertex distances:")
print(sqrt(myhexc))


## @knitr polyex1

library(minqa)
cat("Attempt with quadratic penalty\n")
start <- myhex$par0 # starting parameters (slightly reduced regular hexagon)
lb <- myhex$lb
ub <- myhex$ub
cat("Starting parameters:")
print(start)
sol1 <- bobyqa(start, polyobj1, lower=lb, upper=ub, control=list(iprint=2), penfactor=100)
print(sol1)


## @knitr polyex2

## library(minqa)
cat("Attempt with logarithmic barrier\n")

start <- myhex$par0 # starting parameters (slightly reduced regular hexagon)
lb <- myhex$lb
ub <- myhex$ub
cat("Starting parameters:")
print(start)
sol2 <- bobyqa(start, polyobj2, lower=lb, upper=ub, control=list(iprint=2), penfactor=.001)
print(sol2)

## @knitr polyex2


restart <- start
bestarea <- 0
area <- polyarea(nv, restart)
pf <- 0.01
while (bestarea < area) {
  bestarea <- area
  sol2a <- optim(start, polyobj2a, control=list(trace=2, maxit=10000), penfactor=pf)
  sol2a
  restart <- sol2a$par
  area <- polyarea(nv, restart)
  cat("penfactor = ", pf,"   area = ",area,"\n")
  pf <- pf*0.1
  tmp <- readline("Next cycle")
}

x0 <- myhex$par0

bmeth <- c("nmkb", "hjkb", "bobyqa")

smult <- opm(x0, polyobj, lower=lb, upper=ub, method=bmeth, control=list(trace=1, maxit=10000), penfactor=1e-3)
smult 
