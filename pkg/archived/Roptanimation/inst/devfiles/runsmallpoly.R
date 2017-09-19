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

## @knitr animation1

tmp <- readline("Now try animation")

start <- myhex$par0
pt1 <- PolyTrack$new()

library(minqa)
ub <- c(rep(1,(nv-1)), rep(0.75*pi, (nv-2))) # approx for angles
lb <- c(rep(0, (2*nv-3)))
# sol <- bobyqa(start, polyobj, lower=lb, upper=ub, control=list(rhobeg=1, rhoend=1e-8, iprint=3), penfactor=10)

pf <- 1e-4

while (pf > 1e-12) {
sol <- bobyqa(start, polyobj, lower=lb, upper=ub, control=list(iprint=1), penfactor=pf)
   start <- sol$par
   carea <- polyarea(start)
   cat("area=", carea,"  penfactor = ",pf,"\n")
   pf <- pf / 10
}

tmp <- readline("stop here")
# pt2 <- PolyTrack$new()

# Redo the plots/animation after the optimization
# JN June 3 -- not quite working. Looks like almost.
# tkexamp(pt1$PlotPolys(), 
#    list(i=list('animate', init=1, from=1, to=length(pt1$parms), delay=pt1$Delay*100)),
#    hscale=1.25, vscale=1.25)

tmp <- readline("continue to rest of examples")

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

## @knitr polyexq

start <- myhex$par0 # starting parameters (slightly reduced regular hexagon)
lb <- myhex$lb
ub <- myhex$ub
cat("Starting parameters:")
print(start)

library(minqa)
cat("Attempt with quadratic penalty\n")
sol1 <- bobyqa(start, polyobjq, lower=lb, upper=ub, control=list(iprint=2), penfactor=100)
print(sol1)
cat("area = ",polyarea(sol1$par),"\n")

## @knitr polyexbig

library(optimrx)
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

x0a <- sol2$par
sol2a <- bobyqa(x0a, polyobj, lower=lb, upper=ub, control=list(iprint=2), penfactor=1e-6)
print(sol2a)
cat("Area found=",polyarea(sol2a$par),"\n")

x0b <- sol2a$par
sol2b <- bobyqa(x0b, polyobj, lower=lb, upper=ub, control=list(iprint=2), penfactor=1e-9)
print(sol2b)
cat("Area found=",polyarea(sol2b$par),"\n")

## But a further attempt does very poorly.
x0c <- sol2b$par
sol2c <- bobyqa(x0c, polyobj, lower=lb, upper=ub, control=list(iprint=2), penfactor=1e-12)
print(sol2c)
cat("Area found=",polyarea(sol2c$par),"\n")
```

## @knitr polyex2aa

## library(optimrx)
## cat("Attempt with logarithmic barrier using nmkb and hjkb\n")

## sol2a <- opm(x0, polyobjbig, method=meths, bignum=1e+10)
## print(sol2a)


## @knitr polyex3g

library(Rvmmin)
cat("try to reduce the penalty factor. Rvmmin minimizer on polyobj\n")
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
#  tmp <- readline("Next cycle")
}
cat("Parameters from polyex3g\n")
sol3vpar <- sol3v$par
f <- polyobj(sol3vpar, penfactor=pf)
cat("Objective =", f," area =",attr(f,"area"),"  minslack=",attr(f,"minslack"),"\n")


## @knitr polyex4

x0 <- myhex$par0
bmeth <- c("nmkb", "hjkb", "bobyqa")
library(optimrx)
smult <- opm(x0, polyobj, lower=lb, upper=ub, method=bmeth, control=list(trace=1, maxit=10000), penfactor=1e-3)
print(smult )


## @knitr polyex5

## x0 <- myhex$par0
## library(nloptr)
## cat("Still have to put in nloptr calls\n")


## @knitr polyexuall

library(optimrx)
methset <- c("Rvmmin", "L-BFGS-B", "nlminb")
suall <- opm(x0, polyobju, polygradu, method=methset, control=list(trace=0, kkt=FALSE), penfactor=1e-5)
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

# library(optimrx)
bmeth <- c("bobyqa", "L-BFGS-B", "lbfgsb3", "Rvmmin", "Rtnmin", "Rcgmin", "nlminb", "nmkb", "hjkb", "hjn")
suball <- opm(x0, polyobj, polygrad, lower=lb, upper=ub, method=bmeth, 
        control=list(trace=0, kkt=FALSE), penfactor=1e-5)
# NOTE: Got complex Hessian eigenvalues when trying for KKT tests
suball <- summary(suball, order=value)
print(suball)
resb <- coef(suball)
nmeth <- dim(resb)[1]

## @knitr polyexhjn

library(optimrx)
# repeat of earlier code to ensure we have start and bounds
start <- myhex$par0 # starting parameters (slightly reduced regular hexagon)
lb <- myhex$lb
ub <- myhex$ub
cat("Starting parameters:")
print(start)
x0 <- start

shjnp <- opm(x0, polyobj, polygrad, lower=lb, upper=ub, method="hjn", 
        control=list(trace=0, kkt=FALSE), penfactor=1e-5)
shjnp
tmp <- readline("continue")

shjnp1 <- optimr(x0, polyobj, polygrad, lower=lb, upper=ub, method="hjn", hessian=FALSE,
        control=list(trace=0, kkt=FALSE), penfactor=1e-5)
shjnp1
tmp <- readline("continue")

shjn0p <- hjn(x0, polyobj, lower=lb, upper=ub, bdmsk=NULL, control=list(trace=0), penfactor=1e-5)
shjn0p

## @knitr polyexlbfgs

   newfn <- function(spar, fonly=FALSE,  ...){
      f <- efn(spar, ...)
      if (! fonly) { 
         g <- egr(spar, ...)
         attr(f,"gradient") <- g
      } else { attr(f, "gradient") <- NULL }
      attr(f,"hessian") <- NULL # ?? maybe change later
      f
   }
   efn <- polyobj
   egr <- polygrad # to define

 library(Rtnmin)
 stnb <- tnbc(x0, newfn, lower=lb, upper=ub, trace=TRUE)
 stnb

bmeth <- c("L-BFGS-B", "lbfgsb3", "Rtnmin")
suball <- opm(x0, polyobjp, polygrad, lower=lb, upper=ub, method=bmeth, 
        control=list(trace=1, kkt=FALSE), penfactor=1e-5)
# NOTE: Got complex Hessian eigenvalues when trying for KKT tests
suball <- summary(suball, order=value)
print(suball)
resb <- coef(suball)
nmeth <- dim(resb)[1]

