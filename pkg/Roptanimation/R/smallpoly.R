polysetup <- function(nv, defsize=0.98, qpen=.1){
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

polypar2XY <- function(nv, b) {
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

polyarea<-function(nv, b) {
   # compute area of a polygon defined by radial coordinates
   area <- 0 
   l8 <- nv-3
   for (l in 3:nv){
      q1 <- b[l-2]
      q2 <- b[l-1]
      q3 <- b[l+l8]
      atemp <- q1*q2*sin(q3)
      area <- area + atemp
   }
   area <- area * 0.5
   area
}

polydistXY <- function(nv, XY) {
#   compute point to point distances from XY data
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

# Example code
nv <- 6
cat("Polygon data:\n")
myhex <- polysetup(nv)
print(myhex)
cat("Area:\n")
myhexa <- polyarea(nv, myhex$par0)
print(myhexa)
cat("XY coordinates\n")
myheXY <- polypar2XY(nv, myhex$par0)
print(myheXY)
plot(myheXY$x, myheXY$y, type="l")
cat("Constraints:\n")
myhexc<-polydistXY(nv, myheXY)
print(myhexc)
cat("Vertex distances:")
print(sqrt(myhexc))
