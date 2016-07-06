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
   nv <- (length(pars) + 3)/2
   XY <- polypar2XY(pars)
   dist2 <- polydistXY(XY)
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
 # (negative area) + penfactor*(sum(squared violations))
 nv = (length(x)+3)/2 # number of vertices
 f <-  - polyarea(x) # negative area
 dist2 <- polypardist2(x) # from radial coords, excluding radii (bounded)
 slacks <- 1.0 + epsilon - dist2 # slack vector
 if (any(slacks <= 0)) { f <- bignum } # in case of step into infeasible zone
 else {  f <- f - penfactor*sum(log(slacks)) }
 f
}


## @knitr polygrad

polygrad <- function(x, penfactor=1e-8, epsilon=0) {
 nv <- (length(x)+3)/2
 l8 <- nv - 3 # end of radii params
# epsilon <- 0
 bignum <- 1e+20
 # (negative area) + penfactor*(sum(squared violations))
 nn <- length(x)
 gg <- rep(0, nn)
 dist2 <- polypardist2(x) # from radial coords, excluding radii (bounded)
 slacks <- 1.0 + epsilon - dist2 # slack vector
 if (any(slacks <= 0)) { stop("Infeasible") } 
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


## @knitr polyobju

polyobju <- function(x, penfactor=1e-8, epsilon=0) {
# polyobj with radial parameters constrained by log barrier
# epsilon <- 0
 bignum <- 1e+20
 # (negative area) + penfactor*(sum(squared violations))
 nv = (length(x)+3)/2 # number of vertices
 f <-  - polyarea(x) # negative area
 dist2 <- polypardist2(x) # from radial coords, excluding radii (bounded)
 dist2 <- c(x[1:(nv-1)]^2, dist2) # Note the squared distances used
 slacks <- 1.0 + epsilon - dist2 # slack vector
 if (any(slacks <= 0)) { f <- bignum } # in case of step into infeasible zone
 else {  f <- f - penfactor*sum(log(slacks)) }
 f
}

## @knitr polygradu

polygradu <- function(x, penfactor=1e-8, epsilon=0) {
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


## @knitr polyobjq

polyobjq <- function(x, penfactor=0) {
 # negative area + penfactor*(sum(squared violations))
 nv = (length(x)+3)/2 # number of vertices
 f <-  -polyarea(x) # negative area
 XY <- polypar2XY(x)
 dist2 <- polydistXY(XY)
 viol <- dist2[which(dist2 > 1)] - 1.0
 f <- f + penfactor * sum(viol)
 f
}

## @knitr polyobjbig

polyobjbig <- function(x, bignum=1e10) {
 # Put objective to bignum when constraints violated
 nv = (length(x)+3)/2 # number of vertices
 d2 <- c(x[1:(nv-1)]^2, polypardist2(x)) # distances
 if (any(d2 >=1)) { f <- bignum }
 else { f <-  -polyarea(x) } # negative area
 f
}


## @knitr polyex0

# Example code
nv <- 6
cat("Polygon data:\n")
myhex <- polysetup(nv)
print(myhex)
x0 <- myhex$par0 # initial parameters
cat("Area:\n")
myhexa <- polyarea(x0)
print(myhexa)
cat("XY coordinates\n")
myheXY <- polypar2XY(x0)
print(myheXY)
plot(myheXY$x, myheXY$y, type="l")
cat("Constraints:\n")
myhexc<-polydistXY(myheXY)
print(myhexc)
cat("Vertex distances:")
print(sqrt(myhexc))
cat("check distances with polypar2distXY\n")
try1 <- polypar2distXY(x0)
print(try1)
cat("check distances with polypardist2 augmenting output with parameter squares\n")
try2 <- polypardist2(x0)
try2 <- c(x0[1:(nv-1)]^2, try2)
print(try2)
cat("Max abs difference = ",max(abs(try1-try2)),"\n")



## @knitr polyexq

library(minqa)
cat("Attempt with quadratic penalty\n")
start <- myhex$par0 # starting parameters (slightly reduced regular hexagon)
lb <- myhex$lb
ub <- myhex$ub
cat("Starting parameters:")
print(start)
sol1 <- bobyqa(start, polyobjq, lower=lb, upper=ub, control=list(iprint=2), penfactor=100)
print(sol1)
cat("area = ",polyarea(sol1$par),"\n")

## @knitr polyexbig

library(optimr)
cat("Attempt with setting objective big on violation\n")

x0 <- myhex$par0 # starting parameters (slightly reduced regular hexagon)
cat("Starting parameters:")
print(x0)
meths <- c("Nelder-Mead", "nmkb", "hjkb", "newuoa")
solb <- opm(x0, polyobjbig, method=meths, bignum=1e+10)
print(solb)

## @knitr polyexbigplot

NMpar <- unlist(solb["Nelder-Mead",1:9])
nmkbpar <- unlist(solb["nmkb",1:9])
print(NMpar)
cat("Nelder-Mead area=", polyarea(NMpar))
print(nmkbpar)
cat("nmkb area=", polyarea(nmkbpar))
NMXY <- polypar2XY(NMpar)
nmkbXY <- polypar2XY(nmkbpar)
plot(NMXY$x, NMXY$y, col="red", type="l", xlim=c(-.25,0.85), ylim=c(-.05, 1.05), xlab="x", ylab="y")
points(nmkbXY$x, nmkbXY$y, col="blue", type="l")
title(main="Hexagons from NM (red) and nmkb (blue)")


## @knitr polyex2

library(minqa)
cat("Attempt with logarithmic barrier\n")

x0 <- myhex$par0 # starting parameters (slightly reduced regular hexagon)
lb <- myhex$lb
ub <- myhex$ub
cat("Starting parameters:")
print(x0)
sol2 <- bobyqa(x0, polyobj, lower=lb, upper=ub, control=list(iprint=2), penfactor=1e-3)
print(sol2)
cat("Area found=",polyarea(sol2$par),"\n")

## @knitr polyex2a

## library(optimr)
## cat("Attempt with logarithmic barrier using nmkb and hjkb\n")

## sol2a <- opm(x0, polyobjbig, method=meths, bignum=1e+10)
## print(sol2a)


## @knitr polyex3g

library(Rvmmin)
cat("try to reduce the penalty factor. Rvmmin minimizer on polyobju\n")
restart <- x0
bestarea <- 0
lb <- myhex$lb
ub <- myhex$ub
area <- polyarea(x0)
pf <- 0.01
while (bestarea + 1e-14 < area) {
  bestarea <- area
  sol3v <- Rvmmin(restart, polyobj, polygrad, lower=lb, upper=ub, control=list(trace=2, maxit=1000), penfactor=pf)
  sol3v
  restart <- sol3v$par
  area <- polyarea(restart)
  cat("penfactor = ", pf,"   area = ",area," change=",area-bestarea,"\n")
  pf <- pf*0.1
  tmp <- readline("Next cycle")
}


## @knitr polyex4

x0 <- myhex$par0
bmeth <- c("nmkb", "hjkb", "bobyqa")
library(optimr)
smult <- opm(x0, polyobj, lower=lb, upper=ub, method=bmeth, control=list(trace=1, maxit=10000), penfactor=1e-3)
print(smult )


## @knitr polyex5

## x0 <- myhex$par0
## library(nloptr)
## cat("Still have to put in nloptr calls\n")


## @knitr polyexuall

library(optimr)
suall <- opm(x0, polyobju, polygradu, control=list(all.methods=TRUE, kkt=FALSE), penfactor=1e-5)
# NOTE: Got complex Hessian eigenvalues when trying for KKT tests
suall <- summary(suall, order=value)
print(suall)
resu <- coef(suall)
nmeth <- dim(resu)[1]

## @knitr allplotu

mheXY <- polypar2XY(x0)
plot(mheXY$x, mheXY$y, col='red', type='l', xlim=c(-0.5, 1.05), ylim=c(-0.1, 1.2), xlab='x', ylab='y')
for (ii in 1:nmeth){
   mpar <- resu[ii,]
   XY <- polypar2XY(mpar)
   points(XY$x, XY$y, type='l', col='green')
}
   
## @knitr polyexallb

# library(optimr)
bmeth <- c("bobyqa", "L-BFGS-B", "lbfgsb3", "Rvmmin", "Rtnmin", "Rcgmin", "nlminb", "nmkb", "hjkb")
suball <- opm(x0, polyobj, polygrad, lower=lb, upper=ub, method=bmeth, 
        control=list(kkt=FALSE), penfactor=1e-5)
# NOTE: Got complex Hessian eigenvalues when trying for KKT tests
suball <- summary(suball, order=value)
print(suball)
resb <- coef(suball)
nmeth <- dim(resb)[1]

